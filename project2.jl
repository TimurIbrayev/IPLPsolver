
using Printf
using LinearAlgebra
using MatrixDepot
using DataStructures
using SparseArrays
using AMD

# Authored by: Timur Ibrayev and Gobinda Saha
# Last Modified: 5/4/2020


"""
Data structure for input LP problem
"""
mutable struct IplpProblem
  c::Vector{Float64}
  A::SparseMatrixCSC{Float64} 
  b::Vector{Float64}
  lo::Vector{Float64}
  hi::Vector{Float64}
end

"""
Data structure for storing presolved portion of LP problem
"""
mutable struct IplpPresolved
  xpsol::Vector{Float64}
  cpsol::Vector{Float64}
end

"""
Data structure for storing solution for LP problem using IP solver
"""
mutable struct IplpSolution
  x::Vector{Float64} # the solution vector 
  flag::Bool         # a true/false flag indicating convergence or not
  cs::Vector{Float64} # the objective vector in standard form
  As::SparseMatrixCSC{Float64} # the constraint matrix in standard form
  bs::Vector{Float64} # the right hand side (b) in standard form
  xs::Vector{Float64} # the solution in standard form
  lam::Vector{Float64} # the solution lambda in standard form
  s::Vector{Float64} # the solution s in standard form
end

"""
	free_variable(A,b,c,P_lo,P_hi)

Finds and removes free variables from the given input problem.
"""
function free_variable(A,b,c,P_lo,P_hi)
    idx_lo=findall( x->x.<=-1e308, vec(P_lo))
    idx_hi=findall( x->x.>=1e308, vec(P_hi))
    free_var_idx=collect(intersect!(OrderedSet(idx_hi),OrderedSet(idx_lo)))
    if isempty(free_var_idx)
        @info "No free variable"
        As=A
        cs=c        
    else
        @info "Free variable Found!"
        println("Free variable Number:",length(free_var_idx))
        m,n = size(A)
        fv_len= length(free_var_idx)
        lo=vec(copy(P_lo))
        lo[free_var_idx]=zeros(fv_len)
        P_lo=vec([lo;zeros(fv_len)])
        hi=vec(copy(P_hi))
        P_hi=vec([hi;1e308*ones(fv_len)])
        As=[A -A[:,free_var_idx]]
        cs=vec([c;-c[free_var_idx]])
    end
    return As,b,cs,P_lo,P_hi, free_var_idx
end


"""
	find_singleton_rows(Ak,bk,ck,lb,ub,k,status)

Finds and removes singleton rows from the given input problem.
"""
function find_singleton_rows(Ak,bk,ck,lb,ub,k,status)
    m,n=size(Ak)
    d=[]
    for i=k
        z=collect(findall(Ak[i,:] .!=0))
        append!(d,z)
    end
    xs=vec(bk[k])./vec(sum(Ak[k,:],dims=2))

    col=collect(setdiff(OrderedSet(collect(1:n)),OrderedSet(d)))
    rm_col=collect(setdiff(OrderedSet(collect(1:n)),OrderedSet(col)))
    row=collect(setdiff(OrderedSet(collect(1:m)),OrderedSet(k)))
    cs=ck[rm_col]
    bk=vec(bk)-Ak[:,rm_col]*xs
    Ak=Ak[:,col]
    Ak=Ak[row,:]
    bk=bk[row]
    ck=ck[col]
    if any((xs .>=vec(lb[rm_col])) .!= (xs .<=vec(ub[rm_col])))
        @info "Infeasible solution from singleton Row"
        status= :Infeasible
    end
    lb=lb[col]
    ub=ub[col]

    return Ak,bk,ck,lb,ub,xs,cs,status        
end

"""
	remove_zero_rows(A,b,tol=1e-12)

Finds and removes zero rows from the given input problem.
"""
function remove_zero_rows(A,b,tol=1e-12)
    Ai=copy(abs.(A))
    v =sum(Ai,dims=2)
    idx  =findall( x->x.>tol, vec(v))
    idx_b=findall( x->x.<tol, vec(v))
    Ar=A[idx,:]
    br=b[idx]
    status = :feasible
    if !isempty(idx_b)
        println("Number of zero rows:",length(idx_b))
        if any(b[idx_b] .!= 0)
            status = :Infeasible
            println("There is a zero row in A_eq with a nonzero corresponding entry in b_eq.")
            @info "The problem is infeasible."
        end
    end
    return Ar,br,status
end

"""
	remove_redundancy_sparse(A,rhs)

Finds and removes dependens rows (if any).
"""
function remove_redundancy_sparse(A,rhs)
    tolapiv = 1e-8
    tolprimal = 1e-8
    status = :feasible
    ## Removing zero rows
    A,rhs,status=remove_zero_rows(A,rhs)
    
    if status != :feasible
        return A, rhs, status
    end
    
    m,n=size(A)
    @show m,n,rank(A)
    if rank(A) != m
        @info "Matrix is not full rank after removal of zero rows"
    else
        @info "Matrix is Full rank after removal of zero rows. No need to check for dependent rows!"
        return A, rhs, status
    end
    
    A_orig=A
    
    try    
        # Imp:1 ----- Checking and removing dependent rows -----
        v=collect(1:m)
        b=copy(v)
        k=collect(m+1:m+n)
        d=[]

        A=[Matrix{Float64}(I, m, m) A] 
        e=zeros(m)
        println("Finding Dependent rows, might take a while ---------")
        for i=b
            B=A[:,b]
            e[i]=1
            if i>1
                e[i-1]=0
            end
            pi=B'\e
            js = collect(setdiff(OrderedSet(k),OrderedSet(b)))
            c=abs.(A[:,js]'*pi)
            idx=findall(x->x.>tolapiv, vec(c))
            if !isempty(idx)
                j=js[idx[1]]
                b[i]=j
            else
                push!(d,i)
            end        
        end

        keep=OrderedSet(collect(1:m))
        k=collect(setdiff(keep,OrderedSet(d)))
        println("Alg1: Number of dependent rows:",length(d))
        return A_orig[k, :], rhs[k], status
    
    catch
    
        ## Imp:1 ----- Checking and removing dependent rows B-----by cholesky
        skip_tol = 1.0e-10
        A=copy(A_orig)
        M = A*A' ## M symetric matix
        m = size(M, 1)
        zero_out = ones(m)

        ## ----- Modified Cholesky Factorization with Skipping small pivot columns       
        for i = 1:m-1
            if M[i,i] <= skip_tol
                Mii = 1.0e64
                M[i, i] = Mii
                M[i+1:m, i] .= 0.0
                zero_out[i] = 0.0
            else
                Mii = sqrt(M[i,i])
                M[i, i] = Mii
                M[i+1:m, i] ./= Mii
                M[i+1:m, i+1:m] -= M[i+1:m, i]*M[i+1:m, i]'
            end
        end
        if M[m, m] <= skip_tol
            M[m,m] = 1.0e64
            zero_out[m] = 0.0
        else
            M[m,m] = sqrt(M[m,m])
        end
        k=findall( x->x.==1, vec(zero_out))
        println("Alg2: Number of dependent rows:",(m-length(k)))
        return A_orig[k, :], rhs[k], status
    end
end

"""
	preprocess_LP(A,b,c,P_lo,P_hi,status,data_status,tol=1e-13)

Presolves the given input problem.
"""
function preprocess_LP(A,b,c,P_lo,P_hi,status,data_status,tol=1e-13)
    # array contains solved values of x after preprocessing
    xsol=vec([0.])
    csol=vec([0.])
    ##-----------(a) delete Fixed variable ---------------##
    rm_idx=findall(vec(P_lo).==vec(P_hi))
    kp_idx=findall(vec(P_lo).!=vec(P_hi))
    lo=vec(copy(P_lo))
    hi=vec(copy(P_hi))
    if !isempty(rm_idx)
        @info "Fixed variable detected"
        println("Number of fixed variables:",length(rm_idx))
        As=A[:,kp_idx]
        cs=c[kp_idx]
        xfix=lo[rm_idx]
        bs=vec(b)-A[:,rm_idx]*sparsevec(xfix)
        P_lo=lo[kp_idx]
        P_hi=hi[kp_idx]
        cfix=c[rm_idx]
        xsol=[xsol;xfix]
        csol=[csol;cfix]
        data_status =:changed
    else
        As=A
        cs=c
        bs=b
        P_lo=lo
        P_hi=hi
    end
    
    ##-----------(b) delete Zero Columns ---------------##
    m,n=size(As)
    zrcol_idx=findall(vec(maximum(abs.(As),dims=1)) .< tol)
    keep=collect(setdiff(OrderedSet(collect(1:n)),OrderedSet(zrcol_idx)))
    if !isempty(zrcol_idx)
        @info "Zero Column detected"
        println("Number of zero columns:",length(zrcol_idx))
        if any((cs[zrcol_idx] .<0) .== (P_hi[zrcol_idx] .>=(1e308)))
            @info "Problem Unbounded below"
            status=:Unbounded
            return As,bs,cs,P_lo,P_hi,  IplpPresolved(xsol,csol), status, data_status,rm_idx,kp_idx,zrcol_idx,keep
        end
        
        xzrcol=zeros(length(zrcol_idx)) + (cs[zrcol_idx] .<0).*P_hi[zrcol_idx] + (cs[zrcol_idx] .>0).*P_lo[zrcol_idx]
        As=As[:,keep]
        ct=copy(cs)
        cs=ct[keep]
        czrcol=ct[zrcol_idx]
        P_lo=P_lo[keep]
        P_hi=P_hi[keep]
        xsol=[xsol;xzrcol]
        csol=[csol;czrcol]
        data_status =:changed

    end
    return As,bs,cs,P_lo,P_hi, IplpPresolved(xsol,csol), status, data_status,rm_idx,kp_idx,zrcol_idx,keep
end

"""
	convert_to_standard(P::IplpProblem)

Converts preprocessed input problem into standard form.
"""
function convert_to_standard(P::IplpProblem)
    println("Initial Constraint Matrix size-")
    @show m,n=size(P.A)
    presolver=IplpPresolved(vec([0.]),vec([0.]))
    obj_corr=0.
    idx_ignore=[] # if remains empty no dat corrected 
    ## check feasibility, remove zero rows and dependent rows ###
    status=:feasible
    data_status=:unchanged
    ## -----------------Step-1: Checking for Bound Infeasibility ---------------##
    if any(P.lo .> P.hi)
        @info "Lowerbound exceeds upper bound"
        status=:Infeasible
        return IplpProblem(vec(P.c), P.A, vec(P.b), vec(P.lo), vec(P.hi)),status, data_status, presolver,obj_corr,idx_ignore
    end
    if sum(P.c)==0
        @info "Minimizing constant objective function f(x)=0 ! Constraints does not mean anything."
        status=:Infeasible
        return IplpProblem(vec(P.c), P.A, vec(P.b), vec(P.lo), vec(P.hi)),status, data_status, presolver,obj_corr,idx_ignore
    end
    ## -----------------Step-2: (a)Remove Zero Rows and (b)Dependent Rows ---------------##
    if issparse(P.A)
        As,bs,status=remove_redundancy_sparse(P.A,P.b)
        if status != :feasible
            @info "Trivial Infeasiblity of the LP problem detected"
            return IplpProblem(vec(P.c), As, vec(bs), vec(P.lo), vec(P.hi)),status, data_status, presolver,obj_corr,idx_ignore
        end
    else
        @info "A is not a sparse matrix. Do dense matrix removal."
        status=:sparsityissue
        return IplpProblem(vec(P.c), P.A, vec(P.b), vec(P.lo), vec(P.hi)),status, data_status, presolver,obj_corr,idx_ignore
    end
    println("Constraint Matrix size after row removal-")
    @show m,n=size(As)
    
    ## -----------------Step-3: Remove (c)fixed variable (d)zero column (e)Singleton rows ---------------##
    cs=vec(copy(P.c))
    P_lo=copy(P.lo)
    P_hi=copy(P.hi)
    As,bs,cs,P_lo,P_hi,presolver,status,data_status,rm_idx,kp_idx,zrcol_idx,keep=preprocess_LP(As,bs,cs,P_lo,P_hi,status,data_status)
    println(status)
    println("Constraint Matrix size after preprocessing-")
    @show m,n=size(As)
    idx_fv=[]
    
    ## -----------------Step-4: Checking for standard form ---------------##
    if status ==:feasible
        if all(P_lo .==0) && all(P_hi .==1.0e308) 
            @info "The original problem is already in standard form!" 
            standard=:given
        else
            @info "Converting original problem to standard form!" 
            standard=:converted
            ## -----------------Step-5: Checking for Free variable ---------------##
            As,bs,cs,P_lo,P_hi,idx_fv=free_variable(As,bs,cs,P_lo,P_hi)
            println("Constraint Matrix size before adding slacks-")
            @show m,n=size(As)
            bs = bs-As*P_lo
            h  = P_hi-P_lo
            obj_corr=dot(cs,P_lo)
            data_status=:changed

            if !all(vec(h) .>=1.0e308)
                @info "Std Conversion includes agugmenting linear constraints!" 
                idx=findall( x->x.<1e308, vec(h))
                k=length(idx)
                println("Additional number of constraints for slacks: ",k)
                B1=zeros(m,k)
                B2=Matrix{Float64}(I, k, k)
                B3=zeros(k,n)
                B3[:,idx]=B2
                As=[As B1;
                    B3 B2]
                bs=[bs;h[idx]]
                cs=[cs;zeros(k,1)]
                println("Constraint Matrix size after adding slacks-")
            end
        end
        @show m,n=size(As)
        if rank(As) != m
            @info "Again Searching for zero/dependent rows!" 
            As,bs,status=remove_redundancy_sparse(As,bs)
        end
    end
    return IplpProblem(vec(cs), As, vec(bs), vec(P_lo), vec(P_hi)),status, data_status, 
						presolver, obj_corr, [rm_idx,kp_idx,zrcol_idx,keep,idx_fv]
end

"""
	get_original_x(x,idx_correct,std_problem,presolver)

Post processing function used to get original x from standard form solution.
"""
function get_original_x(x,idx_correct,std_problem,presolver)
    println("Postsolving to get back original x ---")
    rm_idx,kp_idx,zrcol_idx,keep,free_var_idx = idx_correct
    #@show size(free_var_idx),size(zrcol_idx),size(rm_idx)

    # Correction for standardization
    xr=x[1:length(std_problem.lo)]
    xr=xr.+std_problem.lo
    # Correction for free variable
    if !isempty(free_var_idx)
        xt=xr[free_var_idx].-xr[(length(xr)-length(free_var_idx))+1:end]
        xr[free_var_idx]=xt
        xr=xr[1:(length(xr)-length(free_var_idx))]
    end
    # Correction for zero cloumn
    if !isempty(zrcol_idx)
        lz=length(zrcol_idx)
        lxr=length(xr)
        xz=zeros((lz+lxr))
        xz[keep]=xr
        xz[zrcol_idx]=presolver.xpsol[length(presolver.xpsol)-lz+1:end]
    else
        xz=xr
    end
    # Correction for fixed variable
    if !isempty(rm_idx)
        lf=length(rm_idx)
        lxz=length(xz)
        xf=zeros((lf+lxz))
        xf[kp_idx]=xz
        xf[rm_idx]=presolver.xpsol[2:lf+1]
    else
        xf=xz
    end
    #@assert length(xf) == length(problem.c) "size mismatch!"
    #obj_print=@sprintf("%.10E", dot(vec(xf),vec(problem.c)))
    #println("Orignal  Objective function value :  ",obj_print)
    return xf
end

"""
	get_standard_xs_cs(x,c,std_problem,presolver)

Function converts x and c from standard form to the input form.
"""
function get_standard_xs_cs(x,c,std_problem,presolver)
    # Correction for standardization
    xr=x[1:length(std_problem.lo)]
    xr=xr.+std_problem.lo
    xs=[xr;presolver.xpsol[2:end]]
    cs=[c[1:length(std_problem.lo)];presolver.cpsol[2:end]]
    obj_print=@sprintf("%.10E", dot(vec(xs),vec(cs)))
    #println("Standard  Objective function value:  ",obj_print)
    return xs,cs
end

""" 
	starting_point(c, A, b)

Starting points function used to get initial starting points.
"""
function starting_point(c, A, b)
    AAt=A*A'
    xi=AAt\b
    xb=A'*xi #
    lam_i=A*c
    lam_b=AAt\lam_i #
    sb=c-A'*lam_b #
    
    dx=-1.5*minimum(xb)
    dx=max(0,dx)
    
    ds=-1.5*minimum(sb)
    ds=max(0,ds)  
    
    xh=xb.+dx
    sh=sb.+ds
    
    ## This might be needed to stop certain infeasibility
    if all(xh .==0)
        xh .+=eps()
    end
    if all(sh .==0)
        sh .+=eps()
    end
    
    xs=xh'*sh
    dxh=(0.5*xs)./sum(sh)
    dsh=(0.5*xs)./sum(xh)
    
    x=xh.+dxh
    s=sh.+dsh
    lam=lam_b
    return x',s',lam'
end

"""
	min_negative_ratio(num_array, den_array)

Util function used to compute minimum negative ratio.
"""
function min_negative_ratio(num_array, den_array)
    @assert size(num_array) == size(den_array) "Computation of minimum negative ratio is allowed only on vectors of the same length!"
    min = Inf
    for (n, d) in zip(num_array, den_array)
        if d < 0
            neg_ratio = -1.0*(n/d)
            if neg_ratio < min
                min = neg_ratio
            end
        end
    end
    return min
end

"""
	step_normal_eq_form_beta(As, D2, L, p, zero_out, Sinv, rc, rb, rxs)

Computes single step using normal equations form and factorization of A.
"""
function step_normal_eq_form_beta(As, D2, L, p, zero_out, Sinv, rc, rb, rxs)
    RHS_l = -rb - As*D2*rc + As*Sinv*rxs
        
    dl = L'\(L\(RHS_l[p]))
    dl = dl.*zero_out
    dl = dl[invperm(p)]
    ds = -rc - As'*dl
    dx = -Sinv*rxs - D2*ds
    return dx, dl, ds
end

"""
	mod_chol_1_juliainbuilt(input_M)

First Modified Cholesky factorization function.

This version utilizes Julia in-built Cholesky
with addition of perturbing diagonal elements
if the input matrix M is not exactly pos. def.

This function is used for small and well defined problems.
"""
function mod_chol_1_juliainbuilt(input_M)
    m = size(input_M, 1)
    M = copy(input_M)
    beta = 1.0e-3
    minDiag = minimum(diag(M))
    if minDiag > 0
        tau = 0.0
    else
        tau = -minDiag + beta
    end
    for k = 0:1000
        M += spdiagm(0 => vec(ones(m))*tau)
        if isposdef(M)
            C = 0
            try
                C = cholesky(M)
            catch
                println("\t    Mod.Cholesky 1 failed! Switched to Mod.Cholesky 2!")
                return 0, 0, 2
            end
            L = sparse(C.L)
            p = C.p
            return L, p, 1
        else
            tau = max(2*tau, beta)
        end
    end
    println("\t    Mod.Cholesky 1 failed! Switched to Mod.Cholesky 2!")
    return 0, 0, 2
end

"""
	mod_chol_3_withskip(input_M)

Third Modified Cholesky factorization function.

This version utilizes modified Cholesky which
implements skipping procedure for small pivots
(usually occuring towards the end of convergence)
by replacing the small pivot with large value
and the corresponding column elements by zeros.

This function is used for problems involving 
matrices which are not pos. def.
"""
function mod_chol_3_withskip(input_M)
    skip_tol = 1.0e-8
    m = size(input_M, 1)
    M = copy(input_M)
    ## ---------- Using Approximate Minimim Degree (AMD) pkg
    p = amd(M)
    M = M[p, p]
    
    ## ---------- Track which pivot has been skipped
    # The corresponding solutions should be set to 0
    zero_out = ones(m)
    
    ## ----- Modified Cholesky Factorization with Skipping small pivot columns       
    for i = 1:m-1
        if M[i,i] <= skip_tol
            Mii = 1.0e64
            M[i, i] = Mii
            M[i+1:m, i] .= 0.0
            zero_out[i] = 0.0
        else
            Mii = sqrt(M[i,i])
            M[i, i] = Mii
            M[i+1:m, i] ./= Mii
            M[i+1:m, i+1:m] -= M[i+1:m, i]*M[i+1:m, i]'
        end
    end
    if M[m, m] <= skip_tol
        M[m,m] = 1.0e64
        zero_out[m] = 0.0
    else
        M[m,m] = sqrt(M[m,m])
    end
    L = sparse(LowerTriangular(M))
    return  L, p, zero_out
end

"""
	adaptive_cholesky(input_M, modchol_mode)

Top-level function to select which Cholesky factorization to apply.
"""
function adaptive_cholesky(input_M, modchol_mode)
    zero_out = ones(size(input_M, 1))
    if modchol_mode == 1
        L, p, modchol_mode = mod_chol_1_juliainbuilt(input_M)
    end
    if modchol_mode == 2
        L, p, zero_out = mod_chol_3_withskip(input_M)
    end
    return L, p, zero_out, modchol_mode
end

"""
	MPCsolver_beta(cs, As, bs, starting_points, status;
    		tol=1.0e-8, maxiter=100, modified_choleksy=1, quiet=true)

Mehrotra Predictor-Corrector based IP solver.
"""
function MPCsolver_beta(cs, As, bs, starting_points, status;
    tol=1.0e-8, maxiter=100, modified_choleksy=1, quiet=true)
    # first get the sizes
    m, n = size(As) # m rows * n columns
    # check to see if the matrix is wide
    @assert m < n "Input constraint matrix in standard form As should be wide!"
        
    # fetch starting points
    x = starting_points[1:n]
    l = starting_points[n+1:n+m]
    s = starting_points[n+m+1:2n+m]   
    flag = false
    norm_init=norm([x;s],1)
    
    # select initial modified Cholesky mode (defaul = 1)
    @assert (modified_choleksy == 1) | (modified_choleksy == 2) "Currently this script has only 2 types of modified Cholesky!"
    modchol_mode = modified_choleksy
    
    # DEBUG:
    histx = zeros(n, maxiter)
    objective_value=zeros(0)
    suggest_tol =1.
    
    # Display progress, if enabled
    if !quiet
        @printf("  %6s\t%9s\t%9s\t%9s\t%9s\n", 
            "iter", "duality_measure", "norm_residual", "factorization_residual","variable_norm" );
    end
    
    for iter = 1:maxiter
        ## --------------- Cholesky speed-up -----------------##
        # Since x and s do not change within single iteration
        # and also corresponding matrix Dsquared (D2),
        # Cholesky factorization can be computed only once
        # and then re-used in both affine and centered step!
        D2 = spdiagm(0 => vec(x./s))
        M = As*D2*As'
        L, p, zero_out, modchol_mode = adaptive_cholesky(M, modchol_mode)
        fact_res = norm(L*L' - M[p,p])
        Sinv = spdiagm(0 => vec(1.0./s))
        

        ## ----------------- Predictor step -------------------## 
        # Compute F
        rc = As'*l + s - cs
        rb = As*x - bs
        rxs_aff = x.*s
        # Solve for affine-scaling direction
        dx_aff, dl_aff, ds_aff = step_normal_eq_form_beta(As, D2, L, p, zero_out, Sinv, rc, rb, rxs_aff)  
        

        ## ------------------ Centering step -------------------##
        # Compute max allowable steplengths along the affine-scaling direction
        alpha_aff_prim = min(1.0, min_negative_ratio(x, dx_aff))
        alpha_aff_dual = min(1.0, min_negative_ratio(s, ds_aff))
        # Compute duality measures: running average and due to affine-scaling direction
        mu = x'*s/n
        mu_aff = (x + alpha_aff_prim*dx_aff)'*(s + alpha_aff_dual*ds_aff)/n
        # Compute centering parameter
        sigma  = (mu_aff/mu)^3

        
        ## ------------------- Corrector step -------------------##
        # (to be precise, center-corrector step together)
        # Compute F_corrector
        rxs_cor = rxs_aff + dx_aff.*ds_aff .- sigma*mu
        # Solve aggregated system
        dx, dl, ds = step_normal_eq_form_beta(As, D2, L, p, zero_out, Sinv, rc, rb, rxs_cor)

        
        ## ---------------------- Update step --------------------##
        # Compute max allowed steplength to not violate positivity constraints (i.e. (x,s)>0)
        alpha_max_prim = min_negative_ratio(x, dx)
        alpha_max_dual = min_negative_ratio(s, ds)
        # Decide steplength scaling parameter
        nu = 0.9 + 0.1*(iter-1)/maxiter #heuristic from book, which makes nu -> 1 as iter -> maxiter
        alpha_prim = min(1.0, nu*alpha_max_prim)
        alpha_dual = min(1.0, nu*alpha_max_dual)
        # Update primary and dual variables
        x = x + alpha_prim*dx
        l = l + alpha_dual*dl
        s = s + alpha_dual*ds
        
        
        ## ------------------------ Monitor step -------------------##
        # store the objective value per iterations 
        obj_val = vec(cs)'*x
        push!(objective_value, obj_val)
        # Check convergence control parameters
        duality_measure = (x'*s)/n
        norm_residual   = norm([As'*l + s - cs; As*x - bs; x.*s])/norm([bs;cs])
        variable_norm   =norm([x;s],1)
        
        # store x values    
        histx[:, iter] = x      
        # Display progress, if enabled
        if !quiet
            @printf("  %6i\t%9.2e\t%9.2e\t%9.2e\t%9.2e\n", 
                iter, duality_measure, norm_residual, fact_res, variable_norm)
        end       
        
                
        if ((duality_measure + norm_residual) < suggest_tol)
            suggest_tol=(duality_measure + norm_residual)
        end

        # check if the solution is found
        if all(duality_measure <= tol) & all(norm_residual <= tol)
            flag = true
            histx = histx[:, 1:iter]
            break
            
        elseif iter == maxiter
            # Check Infeasibility of the given LP problem
            if (variable_norm > 1e3*norm_init)
                if (duality_measure > tol) | (variable_norm > 1e9*norm_init)
                    @warn "This LP problem is Infeasible."
                    status =:Infeasible
                else
                    @warn "Feasible problem.Failed to converge for given tolerance within maxiter."
                    x=@sprintf("%.2E",suggest_tol )
                    println("Please try tolerance above: ",x)                        
                end
                    
            else
                @warn "Feasible problem.Failed to converge for given tolerance within maxiter."
                x=@sprintf("%.2E",suggest_tol )
                println("Please try tolerance above: ",x) 
            end
            break
        elseif (variable_norm > 1e6*norm_init) & (duality_measure > tol)
            ## greedy early Infeasibility detection 
            @warn "This LP problem is Infeasible."
            status =:Infeasible
            break            
        end

    end
    return [flag, cs, As, bs, x, l, s], status, histx, objective_value
end 

"""
soln = iplp(Problem,tol) solves the linear program:

   minimize c'*x where Ax = b and lo <= x <= hi

where the variables are stored in the following struct:

   Problem.A
   Problem.c
   Problem.b   
   Problem.lo
   Problem.hi

and the IplpSolution contains fields

  [x,flag,cs,As,bs,xs,lam,s]

which are interpreted as   
a flag indicating whether or not the
solution succeeded (flag = true => success and flag = false => failure),

along with the solution for the problem converted to standard form (xs):

  minimize cs'*xs where As*xs = bs and 0 <= xs

and the associated Lagrange multipliers (lam, s).

This solves the problem up to 
the duality measure (xs'*s)/n <= tol and the normalized residual
norm([As'*lam + s - cs; As*xs - bs; xs.*s])/norm([bs;cs]) <= tol
and fails if this takes more than maxit iterations.
"""
function iplp(Problem::IplpProblem, tol::Float64; maxit::Int64=100, quiet::Bool=false)
	@assert tol > 0.0 "Input tolerance should be positive!"
	@assert maxit > 0 "Input maxit should be greater than 0!"
	# Preprocess the given input problem
	my_problem = IplpProblem(vec(Problem.c), Problem.A, vec(Problem.b), vec(Problem.lo), vec(Problem.hi))
	println("-------------- Presolving --------------")
	std_problem, problem_status, data, presolver, obj_corr, idx_correct = convert_to_standard(my_problem)
	@show problem_status, data, obj_corr

	# if problem is not trivially infeasible...
	if problem_status==:feasible
		## Get the starting points 
	    println("--------- Calculating starting point ----------")
		x_init,s_init,l_init=starting_point(std_problem.c, std_problem.A, std_problem.b)
		
		## Solve the IpLP with MPCsolver
		println("--------- Solving using MPC algorithm ---------" )
		solution,problem_status,histx,objective_value=MPCsolver_beta(std_problem.c, std_problem.A, std_problem.b, 
		    [x_init l_init s_init], problem_status, tol=tol, maxiter=maxit, quiet=quiet) 

		## Postsolving 
		if solution[1]==true  # if converged 
		    if data==:changed
		        println("-------------- Postsolving --------------")
		        objective_value_final = objective_value .+ dot(presolver.xpsol,presolver.cpsol).+ obj_corr
		        x= get_original_x(solution[5],idx_correct,std_problem,presolver)
		        xs,cs=get_standard_xs_cs(solution[5],solution[2],std_problem,presolver)
		    else
		        objective_value_final = objective_value
		        x=solution[5]
		        xs=solution[5]
		        cs=solution[2]
		    end
		    final_solution=IplpSolution(vec(x), solution[1], vec(cs),
		                               solution[3],solution[4],vec(xs),solution[6],solution[7])
			println("\n")		    
			println("--------- Solver Summary ---------" )
		    println("Problem Status     : ",problem_status)
		    println("Solution Converged?: ",final_solution.flag)
		    println("Iteration required : ",length(objective_value_final))
		    xp=@sprintf("%.10E", objective_value_final[end])
		    println("Objective function value:  ",xp)
		else
		    final_solution=IplpSolution(vec([]),solution[1], vec(solution[2]),
		                                solution[3],solution[4],vec([]),vec([]),vec([]))
			println("\n")		    
			println("--------- Solver Summary ---------" )
		    println("Problem Status     : ",problem_status)
		    #@show problem_status
		    #@show final_solution.flag
		end
	else
		final_solution=IplpSolution(vec([]),false, vec(std_problem.c),
		                            std_problem.A,vec(std_problem.b),vec([]),vec([]),vec([]))
		println("\n")		    
		println("--------- Solver Summary ---------" )
	    println("Problem Status     : ",problem_status)
		#@show problem_status
		#@show final_solution.flag
	end

	return final_solution
end

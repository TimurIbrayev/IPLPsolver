This is instructions set to run test of Primal-dual Interior-point LP Solver
using Julia programming language (and its corresponding REPL command line):

# Authored by: Timur Ibrayev and Gobinda Saha
# Last Modified: 5/4/2020

0. Set your current working directory to the folder containing project files
[project2.jl, Manifest.toml, Project.toml, and test.jl], 
by following these steps:
	a. cd("./path/to/project_jl_folder")

1. Install all the required packages, by following these steps:
	a. Open Julia REPL by running "julia" in the terminal
	b. Press "]" to get into Package mode
	c. Run "activate ." to activate "project_jl" environment
	*. Note that REPL prompt should change to (project_jl) to
		reflect the change of environment
	*. Optionally, check required packages by running "status"
	d. Run "instantiate" to install packages required by the project
		in the same state given by the project's Manifest.toml file
	e. Press <<backspace>> or <<CTRL+C>> to exit Package mode
	*. Optinoally, in Julia REPL run "using Pkg; Pkg.status()" and verify that 
		the packages listed are the same as the required packages

2. If your problem IS ALREADY in IplpProblem data structure, 
to run IPLP solver on that problem in Julia REPL follow these steps:

	a. Run "include("./project2.jl")"

	b. Call IplpSolver by running:
		"soln = iplp(problem, tol; maxit, quiet)"

		Note: 
		-problem is required input argument specifying the input problem
				  	which should be in IplpProblem data structure. 
					Typeof: IplpProblem, Default: None, Required

		-tol   	 is required input argument to specify precision tolerance,
				  	which also dictates convergence criteria.
					Typeof: Float64, Default: None, Required

		-maxit   is optional input argument to specify max iterations.
					Typeof: Int64,   Default: 100

		-quiet   is optional input argument to display compute progress.
					Typeof: Bool,    Default: false


3. If your problem IS NOT in IplpProblem data structure, 
to run IPLP solver on some example problem in Julia REPL
follow these steps:

	a. Run "include("./test.jl")"
	b. Run "include("./project2.jl")"

	c. If you want to run the test on some <<NAME>> LP problem 
		from LPNetlib, then specify which problem by running
		"name = "<<NAME>>"" (Name of the problem should be in quotes)
	d. Run "md = mdopen(string("*/", name))"
	e. Run "MatrixDepot.addmetadata!(md.data)"
		Note: this step is required due to issue posted on:
		https://github.com/JuliaMatrices/MatrixDepot.jl/issues/34
	f. Run "problem = convert_matrixdepot(md)"

	h. Finally, call IplpSolver by running:
		"soln = iplp(problem, tol; maxit, quiet)"

		Note: 
		-problem is required input argument specifying the input problem
				  	which should be in IplpProblem data structure. 
					Typeof: IplpProblem, Default: None, Required

		-tol   	 is required input argument to specify precision tolerance,
				  	which also dictates convergence criteria.
					Typeof: Float64, Default: None, Required

		-maxit   is optional input argument to specify max iterations.
					Typeof: Int64,   Default: 100

		-quiet   is optional input argument to display compute progress.
					Typeof: Bool,    Default: false

	*. "test.jl" file contains commented example code with example 
		LP problem ("lp_pilotnov") from LPNetlib for your reference.


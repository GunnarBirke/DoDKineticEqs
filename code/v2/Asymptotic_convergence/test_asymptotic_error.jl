using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("./functions.jl")
include("./vizualize.jl")
include("./functionsRK.jl")
include("./check_parameters.jl")
include("./flux_operators.jl")
include("../IMEX_solve_extern.jl")
include("./analysis_tools.jl")


steps, errors = asymptotic_error(;
                epsilon_refinement = 20, a = 1,
                N = 2^4, Tmax = 3.0, deg = 1, CFL_prefac = 0.5, basis_type = GaussLegendre,
                fluxtype = "altlr", J1_type = "upwind_diss_symm", J1_type_A = "upwind_diss_symm",
                cut_cells = [], fix_eta = true, c = 0.4, include_b_h = true,
                TMM = "ARS3",
                CFL_type = "1/dx",
                                )
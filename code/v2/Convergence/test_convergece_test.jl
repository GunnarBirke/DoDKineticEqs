using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("./functions.jl")
include("./vizualize.jl")
include("./functionsRK.jl")
include("./check_parameters.jl")
include("./flux_operators.jl")
include("../IMEX_solve_extern.jl")
include("./analysis_tools.jl")


cut_cells = []
steps, errors = convergence_test(;
                        eq_type = "transport", epsilon = 0.5, a = 1,
                        exprange = [4,8], Tmax = 3.0, deg = 3, CFL_prefac = 0.5, basis_type = GaussLegendre,
                        fluxtype = "altlr", J1_type = "upwind_diss_symm", J1_type_A = "upwind_diss_symm",
                        cut_cells = cut_cells, fix_eta = true, c = 0.4, include_b_h = true,
                        TMM = SSPRK3)
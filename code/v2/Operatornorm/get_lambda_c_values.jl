using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("./../functions.jl")
include("./../vizualize.jl")
include("./../functionsRK.jl")
include("./../check_parameters.jl")
include("./../flux_operators.jl")
include("../../IMEX_solve_extern.jl")
include("./../analysis_tools.jl")


eq_type = "telegraph"
epsilon = 0.5
a = 1
N = 2^4
deg = 0
CFL_prefac = 0.5 # need to define, but as long fix_eta = ture irrelevant
basis_type = GaussLegendre
J1_type_A = "upwind_diss_symm"
cut_cell_refinement = 51
fix_eta = true
include_b_h = true

flux_tuple = [("altlr" , "upwind_classic"), ("altlr", "upwind_diss_symm"), ("central", "central_symm")]
epsilons = [0.5, 0.25, 0.1, 0.01, 0.001]
deg_max = 5

for eq_type in ["telegraph", "heat"]
    for basis_type in [GaussLegendre, LobattoLegendre]
        for (fluxtype, J1_type) in flux_tuple
            output_matrix = zeros((deg_max+1), 4)
            mkpath("./v2/optimal_lambda_c_values/lambda_c_value_results/" * eq_type * "/$(basis_type)/" * J1_type)
            for (ieps, epsilon) in enumerate(epsilons)
                println("Start: " * eq_type * "/$(basis_type)/" * J1_type * "/epsilon=$(epsilon)")
                for (ideg, deg) in enumerate(0:deg_max)
                    println("deg = $(deg):")
                    @time begin
                        ### get the best norm and best value for c
                        best_norm, best_c = minimize_operatornorm(;
                                                    eq_type = eq_type, epsilon = epsilon, a = a,
                                                    N = N, deg = deg, basis_type = basis_type, 
                                                    fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
                                                    cut_cell_refinement = cut_cell_refinement, fix_eta = fix_eta, include_b_h = include_b_h,
                                                    include_cut_cells = false, do_stabilize = false,
                                                    verbose = false
                                                                        )
                        ### get the respective background operator norm

                        cut_cell_factors, opnorms = op_norm_vs_cut_fac(;
                                                    eq_type = eq_type, epsilon = epsilon, a = a,
                                                    N = N, deg = deg, CFL_prefac = CFL_prefac, basis_type = basis_type,
                                                    fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
                                                    cut_cell_refinement = cut_cell_refinement, fix_eta = fix_eta, c = 1.0, include_b_h = include_b_h,
                                                    include_cut_cells = false, do_stabilize = false,
                                                                        )
                        opnorm_quotient = best_norm/opnorms[1]
                        output_matrix[deg+1, :] = [deg, best_c, best_norm, opnorms[1]]
                    end
                end
                writedlm(joinpath(@__DIR__,"./optimal_lambda_c_values/lambda_c_value_results/" * eq_type * "/$(basis_type)/" * J1_type * "/epsilon=$(epsilon).txt"), output_matrix)
            end
        end
    end
end


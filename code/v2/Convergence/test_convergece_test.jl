using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("../functions.jl")
include("../vizualize.jl")
include("../functionsRK.jl")
include("../check_parameters.jl")
include("../flux_operators.jl")
include("../../IMEX_solve_extern.jl")
include("../analysis_tools.jl")


cut_cells = [0.49, 0.3, 0.1, 10^-2, 10^-7]
plot(xscale = :log10, yscale = :log10)
eq_type = "heat"
solution_output = Matrix{Float64}
for deg in [0, 1, 2]
    steps, errors, solution_output = convergence_test(;
                            eq_type = eq_type, epsilon = 0.5, a = 1,
                            exprange = [4,7], Tmax = 1.0, deg = deg, CFL_prefac = 0.5, basis_type = GaussLegendre,
                            fluxtype = "central", J1_type = "central_symm", J1_type_A = "upwind_diss_symm",
                            cut_cells = cut_cells, fix_eta = true, c = 0.4, include_b_h = true,
                            TMM = SSPRK10_4)

    plot!(steps, errors[:,1], label = "polydeg = $(deg), rho")
    if !(eq_type in ["heat", "transport"])
        plot!(steps, errors[:,2], label = "polydeg = $(deg), j")
    end
end
plot!()
#plot(solution_output[:, 1], solution_output[:, 2], label = "approx")
#plot!(solution_output[:, 1], solution_output[:, 3], label = "exact", linestyle = :dashdotdot)
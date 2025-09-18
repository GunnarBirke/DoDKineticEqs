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
cs = [1.0, 0.55, 0.45]
for (ideg,deg) in enumerate([0, 1, 2])
    steps, errors, solution_output = convergence_test(;
                            eq_type = eq_type, epsilon = 0.5, a = 1,
                            exprange = [4,7], Tmax = 1.0, deg = deg, CFL_prefac = 0.2, basis_type = LobattoLegendre,
                            fluxtype = "altlr", J1_type = "upwind_diss_symm", J1_type_A = "upwind_diss_symm",
                            cut_cells = cut_cells, fix_eta = true, c = cs[ideg], include_b_h = true,
                            TMM = SSPRK3)

    if !(eq_type in ["heat", "transport"])
        plot!(steps[1:Int(length(errors)/2)], errors[1:Int(length(errors)/2)], label = "polydeg = $(deg), rho")
        plot!(steps[1:Int(length(errors)/2)], errors[Int(length(errors)/2)+1:end], label = "polydeg = $(deg), j")
    else
        plot!(steps, errors[:], label = "polydeg = $(deg)")
    end
end
plot!()
#plot(solution_output[:, 1], solution_output[:, 2], label = "approx")
#plot!(solution_output[:, 1], solution_output[:, 3], label = "exact", linestyle = :dashdotdot)
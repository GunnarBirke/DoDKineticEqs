using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("functions.jl")
include("vizualize.jl")
include("functionsRK.jl")
include("IMEX_solve_extern.jl")


#### set up problem
exprange = [2 7]
Tmax = 1.0
deg = 2
TMM = IMEXEuler
epsilon = 0.5
deg_CFL = [0.4 0.3 0.2 0.1]
CFL = deg_CFL[deg+1] * 2*epsilon
basis = GaussLegendre
add_reference = true
global problem
errors = zeros(exprange[2]-exprange[1]+1)
steps = zeros(exprange[2]-exprange[1]+1)
for (iN, N) in enumerate(2 .^(exprange[1]:exprange[2]))
    problem = setup_problem_eq("telsin", -pi, pi, N, Tmax = Tmax, a = 1.0, CFL = CFL, bcs = "periodic", epsilon = epsilon);
    problem = include_cut_cell(problem, 0.00001, 3);

    #### discretize in space
    RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, basis, "Upwind",  do_stabilize = true, fix_eta = true, c = 1, fluxtype = "altlr");
    problem = IMEX_solve_extern(problem, TMM = TMM)
    sol = problem["sol"]
    errors[iN] = problem["error"]
    steps[iN] = problem["cellnumber"]
end

plot(steps, errors, xscale = :log10, yscale = :log10, label = "order $(deg + 1)")
if add_reference == true
    plot!(steps, steps .^-(deg+1), label = "reference order $(deg+1)")
end
plot!(title = "$(basis)_p=$(deg)_eps=$(epsilon)_IMEXEuler", xlabel = "spatial steps", ylabel = "error")

#savefig(plot!(),"./code/dod_telegraph_convergence_results/convergence_$(basis)_p=$(deg)_eps=$(epsilon)_IMEXEuler.pdf")


# plot solution at the largest refinement
#=
sol = problem["sol"]
u_exact = problem["u_exact"]
x_d = problem["x_d"]
full_sol = sol[1:Int(size(sol)[1]/2), :] .+sol[Int(size(sol)[1]/2+1):size(sol)[1], :]*epsilon
vizualize(full_sol, x_d, sol.t, clabel = "full", u_exact = u_exact[1:Int(size(u_exact)[1]/2), :] .+ u_exact[Int(size(u_exact)[1]/2+1):size(u_exact)[1], :]*epsilon)
=#
using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("functions.jl")
include("vizualize.jl")
include("functionsRK.jl")
include("IMEX_solve_extern.jl")


#### set up problem
exprange = [2 8]
Tmax = 1.0
deg = 0
TMM = ImExEuler
epsilon = 0.5
deg_CFL = [0.6 0.25 0.125 0.06 0.06 0.06 0.06]
CFL = deg_CFL[deg+1] * epsilon
basis = GaussLegendre
add_reference = true
include_b_h = false

global problem
errors = zeros(exprange[2]-exprange[1]+1)
steps = zeros(exprange[2]-exprange[1]+1)
for (iN, N) in enumerate(2 .^(exprange[1]:exprange[2]))
    problem = setup_problem_eq("telsin", -pi, pi, N, Tmax = Tmax, a = 1.0, CFL = CFL, bcs = "periodic", epsilon = epsilon);
    #problem = include_cut_cell(problem, 0.00001, 3);

    #### discretize in space
    RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, basis, "Upwind",  do_stabilize = true, fix_eta = true, c = 1, fluxtype = "altlr", include_b_h = include_b_h);
    ex_RHS_mat = problem["ex_RHS_mat"]
    im_RHS_mat = problem["im_RHS_mat"]
    sol, u_exact = TMM(problem, RHS_mat)
    x_d = problem["x_d"]
    v = problem["v"]
    errors[iN] = 0
    for i in 1:problem["cellnumber"]
        basis1 = basis(deg)
        for j in 1:deg+1
            # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
            basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
        end
        errors[iN] += integrate((sol[(deg+1)*(i-1)+1:(deg+1)*i, end] - u_exact[(deg+1)*(i-1)+1:(deg+1)*i, end]).^2, basis1.weights.*(v[i+1]-v[i])/2)
    end
    errors[iN] = sqrt(errors[iN])
    steps[iN] = problem["cellnumber"]
end

plot(steps, errors, xscale = :log10, yscale = :log10, label = "order $(deg + 1)")
if add_reference == true
    plot!(steps, steps .^-(deg+1), label = "reference order $(deg+1)")
end
plot!(title = "$(basis)_p=$(deg)_eps=$(epsilon)_$(TMM)", xlabel = "spatial steps", ylabel = "error")

#savefig(plot!(),"./code/dod_telegraph_convergence_results/convergence_$(basis)_p=$(deg)_eps=$(epsilon)_$(TMM).pdf")
string_bh = ""
if include_b_h == true
    string_b_h = "b_h_on"
else
    string_b_h = "b_h_off"
end
#savefig(plot!(),"./code/telegraph_convergence_results_b_h/without_cut_cells(epsilon=$(epsilon))/p=$(deg)/convergence_$(basis)_p=$(deg)_$(TMM)_$(string_b_h)_CFL=eps*$(deg_CFL[deg+1]).pdf")

# plot solution at the largest refinement
#=
sol = problem["sol"]
u_exact = problem["u_exact"]
x_d = problem["x_d"]
full_sol = sol[1:Int(size(sol)[1]/2), :] .+sol[Int(size(sol)[1]/2+1):size(sol)[1], :]*epsilon
vizualize(full_sol, x_d, sol.t, clabel = "full", u_exact = u_exact[1:Int(size(u_exact)[1]/2), :] .+ u_exact[Int(size(u_exact)[1]/2+1):size(u_exact)[1], :]*epsilon)
=#
plot!()
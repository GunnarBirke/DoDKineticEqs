using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles, LaTeXStrings
include("../functionsRK.jl")
include("../functions.jl")
N = 2^4
Tmax = 3.0
deg = 0
basis = GaussLegendre
TMM = ImEx
TMM_str = "ImExEuler"
nmax = 30
epsilons = (1/2) .^(1:nmax)
CFL_coeff = 0.2
with_cut_cells = false
init_cond = "telsin"

# do not change these following parameters
if with_cut_cells == true
    add_size = 4
else
    add_size = 0
end
vertices = Array(range(-pi, pi, length = N + 1))
h = vertices[2]-vertices[1]
a = 1.0
sol_telegraph = Vector{Matrix{Float64}}()
sol_heat = Vector{Matrix{Float64}}()
errors = zeros(nmax)

for (ieps, eps) in enumerate(epsilons)
    println("Step $(ieps) started: epsilon = $(eps)")
    CFL = CFL_coeff/N^2
    problem = setup_problem_eq(init_cond, -pi, pi, N, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = eps);
    if with_cut_cells == true
        problem = include_cut_cell(problem, 0.001, 3);
        problem = include_cut_cell(problem, 0.1, 7);
        problem = include_cut_cell(problem, 0.0000001, 10);
        problem = include_cut_cell(problem, 0.3, 14);
    end
    # assuming we perform the same timesteps in both calculations
    t_len = length(problem["t_d"])
    RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, basis, "Upwind",  do_stabilize = true, fix_eta = true, c = 0.1, fluxtype = "altlr", include_b_h = true);
    if TMM == ImEx
        push!(sol_telegraph, TMM(problem, RHS_mat, get_RK_tableau(TMM_str))[1])
    else
        push!(sol_telegraph, TMM(problem, RHS_mat)[1])
    end
    problem = setup_problem_eq(init_cond, -pi, pi, N, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = 0);
    if with_cut_cells == true
        problem = include_cut_cell(problem, 0.001, 3);
        problem = include_cut_cell(problem, 0.1, 7);
        problem = include_cut_cell(problem, 0.0000001, 10);
        problem = include_cut_cell(problem, 0.3, 14);
    end
    RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, basis, "Upwind",  do_stabilize = true, fix_eta = true, eq_type = "heat", c = 0.1, fluxtype = "altlr", include_b_h = true);
    if TMM == ImEx
        push!(sol_heat, TMM(problem, RHS_mat, get_RK_tableau(TMM_str))[1])
    else
        push!(sol_heat, TMM(problem, RHS_mat)[1])
    end

    x_d = problem["x_d"]
    v = problem["v"]
    for i in 1:problem["cellnumber"]
        basis1 = basis(deg)
        for j in 1:deg+1
            # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
            basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
        end
        # evaluate at t_len, because thats the last timestep done and the remaining entrys are emtpy
        errors[ieps] += integrate((sol_telegraph[ieps][(deg+1)*(i-1)+1:(deg+1)*i, t_len] - sol_heat[ieps][(deg+1)*(i-1)+1:(deg+1)*i, t_len]).^2, basis1.weights.*(v[i+1]-v[i])/2)
    end
    errors[ieps] = sqrt(errors[ieps])
end

if basis == GaussLegendre
    namestr = "GL"
elseif basis == LobattoLegendre
    namestr = "GLL"
end
if with_cut_cells == true
    cutcellstr = "cut_cells"
else
    cutcellstr = "background"
end

plot(epsilons, errors, xscale = :log10, yscale = :log10, xlabel = L"\epsilon", ylabel = L"\|\|u_h^\epsilon-u_h^0\|\|_M", legend = false, title = "$(cutcellstr), p=$(deg), $(namestr), $(TMM), CFL_coeff=$(CFL_coeff)")
#savefig(plot!(), "./code/heat_equation/Plots/$(TMM)/error_to_heat_eq_$(cutcellstr)_p=$(deg)_$(namestr)_N=$(N).pdf")
#k = 20
#plot(sol_telegraph[k][:, end])
#plot!(sol_heat[k][:, end])

#plot(sol_telegraph[k][1:(deg+1)*2^4, end])
#plot!(sol_heat[k][1:(deg+1)*2^4, end])

#errors_dxsq_dep = errors
#errors_eps_dep = errors
#errors_no_dep = errors
#plot(epsilons, errors_dxsq_dep, xscale = :log10, yscale = :log10, label = L"\cdot/\Delta x^2")
#plot!(epsilons, errors_eps_dep, label = L"\cdot\epsilon/\Delta x")
#plot!(epsilons, errors_no_dep, label = L"\cdot 1/\Delta x")


#plot(sol_telegraph[20][:, length(sol_telegraph[1,1,:])])
#plot!(sol_heat[20][:, length(sol_telegraph[1,1,:])])
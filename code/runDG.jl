using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables
include("functions.jl")
include("vizualize.jl")

#### set up problem
N = 100
Tmax = 0.1
deg = 0
epsilon = 1.01
problem = setup_problem_eq("plateu", 0, 1, N, Tmax = Tmax, a = 1.0, CFL = 0.8, bcs = "periodic", epsilon = epsilon);
#problem = include_cut_cell(problem, 0.13, 3);

#### discretize in space
RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, GaussLegendre, "Upwind",  do_stabilize = true, fix_eta = false, c = 1);

tspan = [0.0,Tmax];
x_d = problem["x_d"];
u0_rho = problem["u0"];
u0_j = zeros(length(u0_rho));
#u0_j = problem["u0"]
u0 = vcat(u0_rho, u0_j);
function A!(du, u, p, t)
    du .= RHS_mat*u
end

tsteps = 100;
t_vec=range(first(tspan), last(tspan), length = tsteps)

prob = ODEProblem(A!, u0, tspan);
sol = solve(prob, Kvaerno3(), saveat=range(first(tspan), stop=last(tspan), length=tsteps));

full_sol = sol[1:(deg+1)*N, :] .+sol[(deg+1)*N+1:(deg+1)*2*N, :]*epsilon


# plots:
# plot rho
plot(camera = (50, 35))
#plot!(t_vec, x_vec, sol[1:N, 1:tsteps],  st=:surface)
# plot g
#plot!(t_vec, x_vec, sol[N+1:2*N, 1:tsteps],  st=:surface)
# plot f
#plot!(t_vec, x_d, full_sol,  st=:surface)

vizualize(full_sol, x_d, t_vec)
#vizualize(sol[1:(deg+1)*N, 1:tsteps], x_d, t_vec)
#vizualize(sol[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps], x_d, t_vec)
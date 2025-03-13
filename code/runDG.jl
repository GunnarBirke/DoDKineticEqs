using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("functions.jl")
include("vizualize.jl")

#### set up problem
N = 10
Tmax = 1.0
deg = 0
epsilon = 0.4
problem = setup_problem_eq("telsin", -pi, pi, N, Tmax = Tmax, a = 1.0, CFL = 1.0, bcs = "periodic", epsilon = epsilon);
#problem = include_cut_cell(problem, 0.13, 3);

#### discretize in space
RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, GaussLegendre, "Upwind",  do_stabilize = true, fix_eta = false, c = 1, fluxtype = "altrl");

tspan = [0.0,Tmax];
x_d = problem["x_d"];
#u0_rho = problem["u0"];
u0 = problem["u0"]
#u0_j = cos.(x_d)
#u0_j = zeros(length(u0_rho));
#u0 = vcat(u0_rho, u0_j);
function A!(du, u, p, t)
    #mul!(du, RHS_mat, u)
    du .= RHS_mat*u
end

tsteps = 100;

# Standard ODE solver
prob = ODEProblem(A!, u0, tspan);
sol = solve(prob, Kvaerno3(), saveat=range(first(tspan), stop=last(tspan), length=tsteps));
#sol = solve(prob, Euler(), dt = 0.001);
#sol = solve(prob, Tsit5(), saveat=range(first(tspan), stop=last(tspan), length=tsteps))


# mit altlr fluss belibt die Lösung auch für kleine epsilon wie die große, starre. Mal untersuchen, wie die einzelnen Komponenten aussehen
#=
# ImeX ODE solver
function A_ex!(du, u, p, t)
    Ns = Int(length(u)/2)
    du[1:Ns] .= RHS_mat[1:Ns, 1:end]*u
    du[Ns+1:2*Ns] .= zeros(Ns)
end
function A_im!(du, u, p, t)
    Ns = Int(length(u)/2)
    du[Ns+1:2*Ns] .= RHS_mat[Ns+1:2*Ns, 1:end]*u
    du[1:Ns] .= zeros(Ns)
end

f = SplitFunction(A_im!, A_ex!)
prob = SplitODEProblem(f,u0, tspan)
sol = solve(prob, SBDF2(), dt = 0.001)
(xsteps, tsteps) = size(sol)
=#

t_vec=range(first(tspan), last(tspan), length = tsteps)
full_sol = sol[1:(deg+1)*N, :] .+sol[(deg+1)*N+1:(deg+1)*2*N, :]*epsilon

# Wave exact solution (with u_0_rho = 1/rsinx, u_0_j=cos(x))
r = -2/(1+sqrt(1-4*epsilon^2))
u_exact = zeros((deg+1)*2*N, tsteps);
u_exact[1:(deg+1)*N, 1:tsteps] = [1/r*exp(r*t)*sin(x) for x in x_d, t in range(0.0, Tmax, length = tsteps)];
u_exact[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps] = [exp(r*t)*cos(x) for x in x_d, t in range(0.0, Tmax, length = tsteps)];


# plots:
# plot rho
plot(camera = (230, 35), xlabel = "t", ylabel = "x")
#plot!(t_vec, x_d, sol[1:N, 1:tsteps],  st=:surface)
# plot g
#plot!(t_vec, x_d, sol[N+1:2*N, 1:tsteps],  st=:surface)
# plot f
#plot!(t_vec, x_d, full_sol,  st=:surface)


#vizualize(full_sol, x_d, t_vec, clabel = "full", u_exact = u_exact[1:(deg+1)*N, 1:tsteps] .+ u_exact[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps])
#vizualize(sol[1:(deg+1)*N, 1:tsteps], x_d, t_vec, clabel = "rho", u_exact = u_exact[1:(deg+1)*N, 1:tsteps])
vizualize(sol[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps], x_d, t_vec, clabel = "g", u_exact = u_exact[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps])


# just exact solution
#vizualize(u_exact[1:(deg+1)*N, 1:tsteps] .+ u_exact[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps]*epsilon, x_d, t_vec, clabel = "full")

#=
t_dest = 2
plot(x_d, u_exact[1:(deg+1)*N, 1] .+ epsilon*u_exact[(deg+1)*N+1:(deg+1)*2*N, t_dest], label = "exact")
plot!(x_d, full_sol[1:(deg+1)*N, t_dest], linestyle = :dash, label = "approx")
=#
#display(RHS_mat)
#writedlm(joinpath(@__DIR__,"./RHS_mat.txt"), round.(RHS_mat, digits = 1))
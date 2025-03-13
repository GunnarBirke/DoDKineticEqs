using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("functions.jl")
include("vizualize.jl")

#### set up problem
N = 200
Tmax = 5.0
deg = 1
epsilon = 0.05
problem = setup_problem_eq("sin", 0, 1, N, Tmax = Tmax, a = 1.0, CFL = 0.03, bcs = "periodic", epsilon = epsilon);
#problem = include_cut_cell(problem, 0.13, 3);

#### discretize in space
RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, LobattoLegendre, "Upwind",  do_stabilize = true, fix_eta = false, c = 1, eq_type = "wave", fluxtype = "full");

tspan = [0.0,Tmax];
x_d = problem["x_d"];
#u0_rho = problem["u0"];
#u0_j = cos.(2*pi*x_d)
#u0 = vcat(u0_rho, u0_j);

u0 = problem["u0"];
function A!(du, u, p, t)
    #mul!(du, RHS_mat, u)
    du .= RHS_mat*u
end

tsteps = 100;

# Standard ODE solver
prob = ODEProblem(A!, u0, tspan);
#sol = solve(prob, Kvaerno3(), saveat=range(first(tspan), stop=last(tspan), length=tsteps));
#sol = solve(prob, Euler(), dt = Tmax/(tsteps-1));
#sol = solve(prob, Tsit5(), saveat=range(first(tspan), stop=last(tspan), length=tsteps))
#sol = solve(prob, CarpenterKennedy2N54(williamson_condition = false), dt = 1/160, save_everystep = true); # solve needs some value here but it will be overwritten by the stepsize_callback

sol, u_exact = SSPRK3(problem, RHS_mat)
#sol, u_exact = expl_Euler(problem, RHS_mat)
tsteps = length(sol[1,:]);


t_vec=range(first(tspan), last(tspan), length = tsteps)
full_sol = sol[1:(deg+1)*N, :] .+sol[(deg+1)*N+1:(deg+1)*2*N, :]*epsilon

# Wave equation exact solution (with u_0_rho = sin(2*pi*x), u0_j = cos(2*pi*x)) (siehe LeVeque, Conservation Laws, Seite 61)
u_wave_exact = zeros((deg+1)*2*N, tsteps);
u_wave_exact[1:(deg+1)*N, 1:tsteps] = [1/2*(sin(2*pi*(x-1/epsilon*t))+epsilon*cos(2*pi*(x-1/epsilon*t))+sin(2*pi*(x+1/epsilon*t))-epsilon*cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in range(0.0, Tmax, length = tsteps)];
u_wave_exact[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps] = [1/2*(1/epsilon*sin(2*pi*(x-1/epsilon*t))+cos(2*pi*(x-1/epsilon*t))-1/epsilon*sin(2*pi*(x+1/epsilon*t))+cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in range(0.0, Tmax, length = tsteps)];

#=
vizualize(full_sol, x_d, t_vec, clabel = "full", u_exact = u_wave_exact[1:(deg+1)*N, 1:tsteps] .+ u_wave_exact[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps])
vizualize(sol[1:(deg+1)*N, 1:tsteps], x_d, t_vec, clabel = "rho", u_exact = u_wave_exact[1:(deg+1)*N, 1:tsteps])
vizualize(sol[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps], x_d, t_vec, clabel = "g", u_exact = u_wave_exact[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps])
=#

#
#i = 25000
i = length(problem["t_d"]) # works just for explicit step number
plot(x_d, u_wave_exact[1:(deg+1)*N, i] .+ u_wave_exact[(deg+1)*N+1:(deg+1)*2*N, i]*epsilon, label = "exact")
#plot!(x_d, u0_rho .+ u0_j*epsilon, linestyle = :dash, label = "init_cond")
plot!(x_d, u0[1:(deg+1)*N] .+ u0[(deg+1)*N+1:(deg+1)*2*N]*epsilon, linestyle = :dash, label = "init_cond")
plot!(x_d, full_sol[:,i], label = "approx", linestyle = :dashdotdot)
#

#writedlm(joinpath(@__DIR__,"./RHS_mat_wave.txt"), round.(RHS_mat, digits = 1))
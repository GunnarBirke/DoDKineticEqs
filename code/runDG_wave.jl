using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("functions.jl")
include("vizualize.jl")
include("functionsRK.jl")
include("IMEX_solve_extern.jl")

#### set up problem
N = 2^4
Tmax = 3.0
deg = 5
epsilon = 1.0
CFL = epsilon*0.5/N
#CFL = 0.25/N
TMM = ImExEuler # Just for IMEX extern solver needed
problem = setup_problem_eq("sin", 0,1, N, Tmax = Tmax, a = 1.0, CFL = CFL, bcs = "periodic", epsilon = epsilon);
#=
problem = include_cut_cell(problem, 0.001, 3);
problem = include_cut_cell(problem, 0.1, 7);
problem = include_cut_cell(problem, 0.0000001, 10);
problem = include_cut_cell(problem, 0.3, 14);
=#

#### discretize in space
RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, GaussLegendre, "Upwind",  do_stabilize = true, fix_eta = true, c = 0.1, fluxtype = "altlr", include_b_h = true, eq_type = "wave");
ex_RHS_mat = problem["ex_RHS_mat"]
im_RHS_mat = problem["im_RHS_mat"]


tspan = [0.0,Tmax];
x_d = problem["x_d"];
u0 = problem["u0"]

#=
# Extern explicit solvers
function A!(du, u, p, t)
    #mul!(du, RHS_mat, u)
    du .= RHS_mat*u
end

tsteps = length(problem["t_d"])

# Standard ODE solver (extern)
prob = ODEProblem(A!, u0, tspan);
#sol = solve(prob, Kvaerno3(), saveat=range(first(tspan), stop=last(tspan), length=tsteps));
#sol = solve(prob, Euler(), dt = 0.001);
#sol = solve(prob, Tsit5(), saveat=range(first(tspan), stop=last(tspan), length=tsteps))

# Telegraph exact solution (with u_0_rho = 1/rsinx, u_0_j=cos(x)) !!!!!just needed, if extern explicit solvers are used!!!!!
r = -2/(1+sqrt(1-4*epsilon^2))
u_exact = zeros((deg+1)*2*N, tsteps);
u_exact[1:(deg+1)*N, 1:tsteps] = [1/r*exp(r*t)*sin(x) for x in x_d, t in range(0.0, Tmax, length = tsteps)];
u_exact[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps] = [exp(r*t)*cos(x) for x in x_d, t in range(0.0, Tmax, length = tsteps)];
=#

#=
# Extern IMEX solver
problem = IMEX_solve_extern(problem, TMM = IMEXEuler)
sol = problem["sol"]
u_exact = problem["u_exact"]
=#

# intern
#sol, u_exact = SSPRK3(problem, RHS_mat)
#sol, u_exact = ImExEuler(problem, RHS_mat)
#sol, u_exact = ARS2(problem, RHS_mat)
#sol, u_exact = ARS3(problem, RHS_mat)

# intern arbitrary
Tableau = get_RK_tableau("ARS3")
#sol, u_exact = ImEx(problem, RHS_mat, Tableau)
sol, u_exact = SSPRK3(problem, RHS_mat)

#t_vec=range(first(tspan), last(tspan), length = tsteps)
full_sol = sol[1:Int(size(sol)[1]/2), :] .+sol[Int(size(sol)[1]/2+1):size(sol)[1], :]*epsilon

# plots:
# plot rho
plot(camera = (230, 35), xlabel = "t", ylabel = "x")
#plot!(problem["t_d"], x_d, sol[1:Int(size(sol)[1]/2), :],  st=:surface)
# plot g
#plot!(problem["t_d"], x_d, sol[Int(size(sol)[1]/2+1):size(sol)[1], :],  st=:surface)
# plot f
#plot!(problem["t_d"], x_d, full_sol,  st=:surface)

########################### animated plot ###########################

# full solution
#vizualize(full_sol, x_d, problem["t_d"], clabel = "full", u_exact = u_exact[1:Int(size(u_exact)[1]/2), :] .+ u_exact[Int(size(u_exact)[1]/2+1):size(u_exact)[1], :]*epsilon)

# 1. component
#vizualize(sol[1:Int(size(sol)[1]/2), :], x_d, problem["t_d"], clabel = "rho", u_exact = u_exact[1:Int(size(u_exact)[1]/2), :])

# 2. component
#vizualize(sol[Int(size(sol)[1]/2)+1:size(sol)[1], :], x_d, problem["t_d"], clabel = "g", u_exact = u_exact[Int(size(u_exact)[1]/2)+1:size(u_exact)[1], :])


########################### plot at a fixed timestep t_dest ###########################
#=
t_dest = size(sol)[2]

#full solution
#plot(x_d, u_exact[1:Int(size(u_exact)[1]/2), t_dest] .+ epsilon*u_exact[Int(size(u_exact)[1]/2+1):size(u_exact)[1], t_dest], label = "exact")
#plot!(x_d, full_sol[:, t_dest], linestyle = :dash, label = "approx")
#
# 1. component
plot(x_d, u_exact[1:Int(size(u_exact)[1]/2), t_dest] , label = "exact")
plot!(x_d, sol[1:Int(size(sol)[1]/2), t_dest], linestyle = :dash, label = "approx")
#
# 2. component
#plot(x_d, u_exact[Int(size(u_exact)[1]/2+1):size(u_exact)[1], t_dest], label = "exact")
#plot!(x_d, sol[Int(size(sol)[1]/2)+1:size(sol)[1], t_dest], linestyle = :dash, label = "approx")
=#



#display(RHS_mat)
#writedlm(joinpath(@__DIR__,"./RHS_mat.txt"), round.(RHS_mat, digits = 1))

# Check eigenvalues for L^2 stability
#=
A = splitRHS["A"]
A_l = size(A)[1]
SBP_A = zeros(A_l*2, A_l*2)
SBP_A[1:A_l, A_l + 1:2*A_l] = A
SBP_A[A_l + 1:2*A_l,1:A_l] = -A

Mh = splitRHS["Mh"]
Dplus = splitRHS["Dplus"]
Dminus = splitRHS["Dminus"]
SBP_B = zeros(A_l*2, A_l*2)
SBP_B[1:A_l, A_l + 1:2*A_l] = Mh*Dminus
SBP_B[A_l + 1:2*A_l,1:A_l] = Mh*Dplus

SBP_B = zeros(A_l*2, A_l*2)
SBP_B[1:A_l, A_l + 1:2*A_l] = Mh*Dminus
SBP_B[A_l + 1:2*A_l,1:A_l] = Mh*Dplus

SBP_C = zeros(A_l*2, A_l*2)
SBP_C[1:A_l, A_l + 1:2*A_l] = -Mh*Dminus
SBP_C[A_l + 1:2*A_l,1:A_l] = -1/epsilon^2*Mh*Dplus 
SBP_C[A_l + 1:2*A_l,A_l + 1:2*A_l] = - 1/epsilon^2*Mh + 1/epsilon*Mh*A
SBP_C = SBP_C


#Mat = SBP_A+SBP_A'å
#Mat = SBP_B+SBP_B'
#Mat = SBP_C + SBP_C'
=#
#
Mat = problem["M"]*RHS_mat + (problem["M"]*RHS_mat)'
#Mat = RHS_mat + RHS_mat'
for ev in eigvals(Mat)
    if real(ev)>10^-12 || abs(imag(ev))>10^-12
        #println(ev)
    end
end
display(eigvals(Mat))
#
#display(RHS_mat)
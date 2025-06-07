using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("../functions.jl")
include("../vizualize.jl")
include("../functionsRK.jl")
include("../IMEX_solve_extern.jl")

#### set up problem
N = 2^4
Tmax = 3.0
deg = 1
epsilon = 0
CFL = 0.1/N
TMM = IMEXEuler # Just for IMEX extern solver needed
problem = setup_problem_eq("telsin", -pi, pi, N, Tmax = Tmax, a = 1.0, CFL = CFL, bcs = "periodic", epsilon = epsilon);
#
problem = include_cut_cell(problem, 0.001, 3);
problem = include_cut_cell(problem, 0.1, 7);
problem = include_cut_cell(problem, 0.0000001, 10);
problem = include_cut_cell(problem, 0.3, 14);
#

#### discretize in space
RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, GaussLegendre, "Upwind",  do_stabilize = true, fix_eta = true, eq_type = "heat", c = 0.1, fluxtype = "altlr", include_b_h = true);


tspan = [0.0,Tmax];
x_d = problem["x_d"];
u0 = problem["u0"]


# intern
#sol, u_exact = SSPRK3(problem, RHS_mat)
sol, u_exact = ImExEuler(problem, RHS_mat)
#sol, u_exact = ARS2(problem, RHS_mat)
#sol, u_exact = ARS3(problem, RHS_mat)

t_dest =  size(sol)[2]

plot(x_d, sol[:, t_dest], label = "approx")
plot!(x_d, u_exact[:, t_dest], linestyle = :dash, label = "u_exact")


Mat = problem["M"]*RHS_mat + (problem["M"]*RHS_mat)'
display(eigvals(Mat))


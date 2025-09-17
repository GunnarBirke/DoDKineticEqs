using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("./functions.jl")
include("./vizualize.jl")
include("./functionsRK.jl")
include("./check_parameters.jl")
include("./flux_operators.jl")
include("./analysis_tools.jl")
include("../IMEX_solve_extern.jl")

##################################
###### set up equation ###########
##################################
# choose initial condition/domain set
#id_set = ("telegraph", "telsin", -pi, pi)
#id_set = ("telegraph_symm", "telsin", -pi, pi)
#id_set = ("telegraph","eq_data", 0, pi)
#id_set = ("wave", "sin", 0, 1)
#id_set = ("wave_symm", "sin", 0, 1)
#id_set = ("transport", "sin", 0, 1) 
id_set = ("heat", "telsin" , -pi, pi)

# equation parameters
epsilon = 0.5
a = 1.0

##################################
###### Numerical parameters ######
##################################
N = 2^4
Tmax = 3.0
deg = 3
CFL = epsilon*0.5/N/(2*deg+1)
fluxtype = "altlr"
J1_type = "upwind_diss_symm"
#fluxtype = "central"
#J1_type = "central_symm"
J1_type_A = "upwind_symm"

##################################
###### set up problem ############
##################################
problem = setup_problem_eq(id_set[1], id_set[2], id_set[3], id_set[4], N, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = epsilon);
#
problem = include_cut_cell(problem, 0.001, 3);
problem = include_cut_cell(problem, 0.4, 7);
problem = include_cut_cell(problem, 0.0000001, 10);
problem = include_cut_cell(problem, 0.1, 14);
#

##################################
###### discretize in space #######
##################################
RHS_mat, problem, SBP_storage = DGsemidiscretization_DoD_telegraph(problem, deg, GaussLegendre,  do_stabilize = true, fix_eta = true,
                                                                 c = 0.5, fluxtype = fluxtype, include_b_h = true, ext_test_func = true,
                                                                 J1_type = J1_type, J1_type_A = J1_type_A);
ex_RHS_mat = problem["ex_RHS_mat"]
im_RHS_mat = problem["im_RHS_mat"]


display(eigvals(problem["M"]*RHS_mat + (problem["M"]*RHS_mat)'))
#vscodedisplay(RHS_mat, "RHS_mat")

##################################
###### discretize in time ########
##################################
tspan = [0.0,Tmax];
x_d = problem["x_d"];
u0 = problem["u0"]
Nx_plot = determine_Nx_plot(problem)

Tableau = get_RK_tableau("ARS3")
#sol, u_exact = ImEx(problem, RHS_mat, Tableau, only_explicit = false)
sol, u_exact = SSPRK3(problem, RHS_mat)
#t_vec=range(first(tspan), last(tspan), length = tsteps)
full_sol = sol[1:(size(sol)[1]÷ 2), :] .+sol[Int(size(sol)[1]÷2+1):size(sol)[1], :]*epsilon
#solimex = sol
solex = sol

##################################
########### plotting #############
##################################
#### animated plot ####
# full solution
#vizualize(full_sol, x_d, problem["t_d"], clabel = "full", u_exact = u_exact[1:Int(size(u_exact)[1]/2), :] .+ u_exact[Int(size(u_exact)[1]/2+1):size(u_exact)[1], :]*epsilon)

# 1. component
#vizualize(sol[1:Nx_plot, :], x_d, problem["t_d"], clabel = "rho", u_exact = u_exact[1:Nx_plot, :], steplim = 500)

# 2. component
#vizualize(sol[Nx_plot+1:size(sol)[1], :], x_d, problem["t_d"], clabel = "g", u_exact = u_exact[Int(size(u_exact)[1]/2)+1:size(u_exact)[1], :])


#### plot at a fixed timestep t_dest #####
#
t_dest = size(sol)[2]#-1000

#full solution
#plot(x_d, u_exact[1:Int(size(u_exact)[1]/2), t_dest] .+ epsilon*u_exact[Int(size(u_exact)[1]/2+1):size(u_exact)[1], t_dest], label = "exact")
#plot!(x_d, full_sol[:, t_dest], linestyle = :dash, label = "approx")
#
# 1. component
plot(x_d, u_exact[1:Nx_plot, t_dest] , label = "exact")
plot!(x_d, sol[1:Nx_plot, t_dest], linestyle = :dash, label = "approx", title = "t = $(problem["t_d"][t_dest])")

# 2. component
if !(id_set[1] in ["transport", "heat"])
    plot!(x_d, u_exact[Int(size(u_exact)[1]/2+1):size(u_exact)[1], t_dest], label = "exact")
    plot!(x_d, sol[Nx_plot+1:size(sol)[1], t_dest], linestyle = :dash, label = "approx")
end
plot!()



##################################
##### some numerical analysis ####
##################################

#=
Dplus = SBP_storage["Dplus"]
Dminus = SBP_storage["Dminus"]
println(zeros(size(Dplus))==round.((Dplus+Dminus'), digits = 12))
display(eigvals((Dplus-Dminus)+(Dplus-Dminus)'))
plot!()
=#


#display(eigvals(dissm + (dissm)'))
#vscodedisplay(dissm, "dissm")
#vscodedisplay(dissp, "dissp")
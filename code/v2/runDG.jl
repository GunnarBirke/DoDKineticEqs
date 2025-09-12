using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("./functions.jl")
include("./vizualize.jl")
include("./functionsRK.jl")
include("./check_parameters.jl")
include("./flux_operators.jl")
include("../IMEX_solve_extern.jl")

#### set up problem
N = 2^4
Tmax = 3.0
deg = 2
epsilon = 0.5
CFL = epsilon*0.5/N
#fluxtype = "altlr"
#J1_type = "upwind_diss_symm"
fluxtype = "central"
J1_type = "central_symm"
J1_type_A = "upwind_diss_symm"
# choose initial condition/domain set
id_set = ("telegraph_symm", "telsin", -pi, pi)
#id_set = ("telegraph","eq_data", 0, pi)
#id_set = ("wave", "sin", 0, 1)
#id_set = ("wave_symm", "sin", 0, 1)
#id_set = ("transport", "sin", 0, 1) 
#id_set = ("heat", "telsin" , -pi, pi)
TMM = ImExEuler # Just for IMEX extern solver needed
problem = setup_problem_eq(id_set[1], id_set[2], id_set[3], id_set[4], N, Tmax = Tmax, a = 1.2, CFL = CFL, bcs = "periodic", epsilon = epsilon);
#
problem = include_cut_cell(problem, 0.001, 3);
problem = include_cut_cell(problem, 0.4, 7);
problem = include_cut_cell(problem, 0.0000001, 10);
problem = include_cut_cell(problem, 0.1, 14);
#
#### discretize in space
RHS_mat, problem, splitRHS, SBP_storage = DGsemidiscretization_DoD_telegraph(problem, deg, GaussLegendre,  do_stabilize = true, fix_eta = true,
                                                                 c = 0.4, fluxtype = fluxtype, include_b_h = true, ext_test_func = true,
                                                                 J1_type = J1_type, J1_type_A = J1_type_A);
#ex_RHS_mat = problem["ex_RHS_mat"]
#im_RHS_mat = problem["im_RHS_mat"]

display(eigvals(problem["M"]*RHS_mat + (problem["M"]*RHS_mat)'))
#vscodedisplay(RHS_mat, "RHS_mat")

tspan = [0.0,Tmax];
x_d = problem["x_d"];
u0 = problem["u0"]
Nx_plot = determine_Nx_plot(problem)

# intern
#sol, u_exact = SSPRK3(problem, RHS_mat)
#sol, u_exact = ImExEuler(problem, RHS_mat)
#sol, u_exact = ARS2(problem, RHS_mat)
#sol, u_exact = ARS3(problem, RHS_mat)

# intern arbitrary
Tableau = get_RK_tableau("ARS3")
#sol, u_exact = ImEx(problem, RHS_mat, Tableau, only_explicit = false)
sol, u_exact = SSPRK3(problem, RHS_mat)
#t_vec=range(first(tspan), last(tspan), length = tsteps)
full_sol = sol[1:(size(sol)[1]÷ 2), :] .+sol[Int(size(sol)[1]÷2+1):size(sol)[1], :]*epsilon

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
#vizualize(sol[1:Nx_plot, :], x_d, problem["t_d"], clabel = "rho", u_exact = u_exact[1:Nx_plot, :], steplim = 500)

# 2. component
#vizualize(sol[Nx_plot+1:size(sol)[1], :], x_d, problem["t_d"], clabel = "g", u_exact = u_exact[Int(size(u_exact)[1]/2)+1:size(u_exact)[1], :])


########################### plot at a fixed timestep t_dest ###########################
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
#



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


#Mat = SBP_A+SBP_A'
#Mat = SBP_B+SBP_B'
#Mat = SBP_C + SBP_C'

=#
#=
Mat = problem["M"]*RHS_mat + (problem["M"]*RHS_mat)'
display(eigvals(Mat))
=#

#display((splitRHS["ex_B"] +  splitRHS["im_B"])[4:7, 4:7])
#display((splitRHS["ex_B_J0"] + splitRHS["im_B_J0"])[4:7, 4:7])
#display((splitRHS["ex_B_J0E1"] + splitRHS["im_B_J0E1"])[4:7, 4:7])
#display((splitRHS["ex_B_J0E2"] + splitRHS["im_B_J0E2"])[4:7, 4:7])
#splitRHS["im_V_J0E1"]
#splitRHS["im_V_J0E2"] = im_V_J0E2

#=
#display(splitRHS["TranspV"]*(splitRHS["ex_B"]+splitRHS["im_B"])*splitRHS["V"])
q1 = 1:18
#q1 = 19:36
q2 = 1:18
#q2 = 19:36

display(splitRHS["VOLTERMS"][q1, q2])
display((-splitRHS["Sc"]*splitRHS["invM"]*(splitRHS["TranspV"]*(splitRHS["ex_B"]+splitRHS["im_B"])*splitRHS["V"] +splitRHS["TranspV"]*(splitRHS["ex_B_J0"]+splitRHS["im_B_J0"])*splitRHS["V"]
    +splitRHS["TranspV"]*(splitRHS["ex_B_J0E1"])*splitRHS["ex_V_J0E1"]+splitRHS["TranspV"]*(splitRHS["im_B_J0E1"])*splitRHS["im_V_J0E1"]
    +splitRHS["TranspV"]*(splitRHS["ex_B_J0E2"])*splitRHS["ex_V_J0E2"]+splitRHS["TranspV"]*(splitRHS["im_B_J0E2"])*splitRHS["im_V_J0E2"]))[q1, q2])
display(RHS_mat[q1, q2])
display(maximum(abs.(RHS_mat)))
=#


#=
#q1 = 1:20
q1 = 21:40
#q2 = 1:20
q2 = 21:40
display(splitRHS["ex_V_J0E1"][q1, q2])
=#
#display(splitRHS["Dminus"])
#display(splitRHS["Dplus"])
#display(problem["x_d"])
#=
println("Dminflux")
display( -splitRHS["Dminusflux"])
println("Dminvol")
display( -splitRHS["Dminusvol"])
println("Dplusflux")
display( -splitRHS["Dplusflux"])
println("Dplusvol")
display( -splitRHS["Dplusvol"])
=#
#display(splitRHS["symm_DM_J1"])
#display(RHS_mat)


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
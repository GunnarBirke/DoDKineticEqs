using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("./../functions.jl")
include("./../vizualize.jl")
include("./../functionsRK.jl")
include("./../check_parameters.jl")
include("./../flux_operators.jl")
include("../../IMEX_solve_extern.jl")
include("./../analysis_tools.jl")

eq_type = "telegraph"
epsilon = 0.5
a = 1
N = 2^4
deg = 0
CFL_prefac = 0.5
basis_type = GaussLegendre
#fluxtype = "altlr"
#J1_type = "upwind_diss_symm"
fluxtype = "central"
J1_type = "central_symm"
J1_type_A = "upwind_diss_symm"
cut_cell_refinement = 51
fix_eta = true
c = 0.4
include_b_h = true



cut_cell_factors, opnorms = op_norm_vs_cut_fac(;
    eq_type = eq_type, epsilon = epsilon, a = a,
    N = N, deg = deg, CFL_prefac = CFL_prefac, basis_type = basis_type,
    fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
    cut_cell_refinement = cut_cell_refinement, fix_eta = fix_eta, c = c, include_b_h = include_b_h,
    include_cut_cells = true, do_stabilize = true,
    )

plot(cut_cell_factors, opnorms, label = "stabilized")
ylim_high = maximum(opnorms)*1.2

cut_cell_factors, opnorms = op_norm_vs_cut_fac(;
    eq_type = eq_type, epsilon = epsilon, a = a,
    N = N, deg = deg, CFL_prefac = CFL_prefac, basis_type = basis_type,
    fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
    cut_cell_refinement = cut_cell_refinement, fix_eta = fix_eta, c = c, include_b_h = include_b_h,
    include_cut_cells = false, do_stabilize = false,
    )

plot!(cut_cell_factors, opnorms, label = "background")
ylim_low = minimum(opnorms)*0.8

cut_cell_factors, opnorms = op_norm_vs_cut_fac(;
    eq_type = eq_type, epsilon = epsilon, a = a,
    N = N, deg = deg, CFL_prefac = CFL_prefac, basis_type = basis_type,
    fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
    cut_cell_refinement = cut_cell_refinement, fix_eta = fix_eta, c = c, include_b_h = include_b_h,
    include_cut_cells = true, do_stabilize = false,
    )

plot!(cut_cell_factors, opnorms, label = "unstabilized")

plot!(ylims = [ylim_low, ylim_high])

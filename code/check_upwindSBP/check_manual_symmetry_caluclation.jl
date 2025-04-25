# just a check that the manual calculations are corret
include("upwindSPB_test.jl")

deg = 4
basis = GaussLegendre
alpha = 0.3
show_matrices = true
digit_tolerance = 12
c = 0.4

IBP, sym, eigenvalues = upwindSBP_test( N = 2^3, 
                                        deg = deg, 
                                        epsilon = 0.5, 
                                        basis = basis, 
                                        c = c, 
                                        digit_tolerance = digit_tolerance, 
                                        alphas = alpha,
                                        print_result = false)                                      


eta = 1-alpha/c
basis_test = basis(deg)
cutcellnodes = 1 .+ alpha .+ basis_test.nodes*alpha
Vol_interpol = interpolation_matrix(cutcellnodes, basis_test)
R = interpolation_matrix(1, basis_test)
L = interpolation_matrix(-1, basis_test)
BJ0 = interpolation_matrix(1+2*alpha, basis_test)
tildeBJ0 = interpolation_matrix((-1-2*alpha/(1-alpha)), basis_test)
display(BJ0)
display(tildeBJ0)

# Symmetry terms
println("--------------------------------------------------------")
println("--------------------------------------------------------")
println("Symmetry terms")
sym_error = (sym-sym')[deg+2:2*deg+2,deg+2:2*deg+2]
analytic_error_term = (-eta*Vol_interpol'*(basis_test.D'*diagm(basis_test.weights)-diagm(basis_test.weights)*basis_test.D)*Vol_interpol)
if show_matrices
    println("Numerical error:")
    display(sym_error)
    println("Analytical error:")
    display(analytic_error_term)
    println("Discrepancy in error terms (rounded up to digit $(digit_tolerance)):")
    display(round.(sym_error - analytic_error_term, digits = digit_tolerance))
end
if round.(sym_error-analytic_error_term, digits = digit_tolerance) == zeros(size(sym_error))
    println("-------Theory matches numerics!-------")
else
    println("-------There is a discrepancy!-------")
end
println("--------------------------------------------------------")
println("--------------------------------------------------------")

# Integration by Parts terms for the generalized setting
println("--------------------------------------------------------")
println("--------------------------------------------------------")
println("IPB terms for the generalized setting")
ibp_error = IBP[(deg+1)+1:2*(deg+1),3*(deg+1)+1:4*(deg+1)]
analytic_error_term_ibp = eta*(R'*tildeBJ0 - BJ0'*L)
if show_matrices
    display(IBP)
    println("Numerical error:")
    display(ibp_error)
    println("Analytical error:")
    display(analytic_error_term_ibp)
    println("Discrepancy in error terms (rounded up to digit $(digit_tolerance)):")
    display(round.(ibp_error - analytic_error_term_ibp, digits = digit_tolerance))
end
if round.(ibp_error - analytic_error_term_ibp, digits = digit_tolerance) == zeros(size(sym_error))
    println("-------Theory matches numerics!-------")
else
    println("-------There is a discrepancy!-------")
end

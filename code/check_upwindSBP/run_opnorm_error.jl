include("upwindSPB_test.jl")

alpha_vec = 0.5 .^range(1, 10, step = 0.1)
IntByParts_opnorm = zeros(length(alpha_vec))
sym_opnorm = zeros(length(alpha_vec))
for (ialpha, alpha) in enumerate(alpha_vec)
    IBP, sym, eigenvalues = upwindSBP_test( N = 2^3, 
                                            deg = 4, 
                                            epsilon = 0.5, 
                                            basis = LobattoLegendre, 
                                            c = 0.4, 
                                            digit_tolerance = 12, 
                                            alphas = alpha,
                                            print_result = false)
    IntByParts_opnorm[ialpha] = opnorm(IBP)
    sym_opnorm[ialpha] = opnorm(sym-sym')
    #
    # just a check that the manual calculations are corret
    if ialpha == Int(length(alpha_vec)/2+0.5)
        deg = 4
        display((sym-sym')[deg+2:2*deg+2,deg+2:2*deg+2])
        eta = 1-alpha/0.4
        basis_test = LobattoLegendre(deg)
        cutcellnodes = 1 .+ alpha .+ basis_test.nodes*alpha
        Vol_interpol = interpolation_matrix(cutcellnodes, basis_test)
        R = interpolation_matrix(1, basis_test)
        L = interpolation_matrix(-1, basis_test)
        display((-eta*Vol_interpol'*(basis_test.D'*diagm(basis_test.weights)-diagm(basis_test.weights)*basis_test.D)*Vol_interpol))
    end
    #
end

plot(alpha_vec, sym_opnorm, label = "symm")
plot!(alpha_vec, IntByParts_opnorm, label = "IntByParts", linestyle = :dashdotdot)



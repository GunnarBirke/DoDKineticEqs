using PolynomialBases, LinearAlgebra

for deg in range(0, 10)
    basis = GaussLegendre(deg)
    alpha = 0.5
    cutcellnodes = 1 .+ 2*alpha .+ basis.nodes*alpha
    Vol_interpol = interpolation_matrix(cutcellnodes, basis)
    println(rank(Vol_interpol) == deg + 1)
end


basis = GaussLegendre(12)
alpha = 0.5
cutcellnodes = 1 .+ 2*alpha .+ basis.nodes*alpha
Vol_interpol = interpolation_matrix(cutcellnodes, basis)
inv(Vol_interpol)
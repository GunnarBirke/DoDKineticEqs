include("upwindSPB_test.jl")

IBP, sym, eigenvalues, problem, splitRHS = upwindSBP_test( N = 2^3, 
                                        deg = 2, 
                                        epsilon = 0.5, 
                                        basis = GaussLegendre, 
                                        c = 0.4,
                                        ext_test_func = false,
                                        digit_tolerance = 12, 
                                        alphas = [0.2],
                                        print_result = true);

#round.(IBP, digits = 12)#+round.(IBP, digits = 12)'
#display(splitRHS["Mh"]*splitRHS["Dplus"])
#display(splitRHS["Dminus"]*splitRHS["Mh"])

IBP;

# check, how the single operators evolve to the result
#=
Mh = splitRHS["Mh"]
display(round.((splitRHS["Dplus"])[1:9, 1:9], digits = 13))
display(round.((splitRHS["Dminus"])[1:9, 1:9], digits = 13))
display(round.((Mh*splitRHS["Dplus"])[1:9, 1:9], digits = 13))
display(round.((splitRHS["Dminus"]'*Mh)[1:9, 1:9], digits = 13))
display(round.((IBP)[1:9, 1:9], digits = 13))
display(round.((IBP+IBP')[1:9, 1:9], digits = 13))
=#

# check, how the background terms add up inside the summation-by-parts structure
#=
xrange = problem["xrange"]
epsilon = problem["epsilon"]
display(( splitRHS["ex_DM"])[1:Int(xrange/2), Int(xrange/2)+1:xrange])
display(( epsilon^2*splitRHS["im_DM"])[Int(xrange/2)+1:xrange, 1:Int(xrange/2)]')
display(( splitRHS["TranspV*ex_B*V"])[1:Int(xrange/2), Int(xrange/2)+1:xrange])
display(( epsilon^2*splitRHS["TranspV*im_B*V"])[Int(xrange/2)+1:xrange, 1:Int(xrange/2)]')
display(( splitRHS["ex_DM"][1:Int(xrange/2), Int(xrange/2)+1:xrange]+epsilon^2*splitRHS["im_DM"][Int(xrange/2)+1:xrange, 1:Int(xrange/2)]'
-splitRHS["TranspV*ex_B*V"][1:Int(xrange/2), Int(xrange/2)+1:xrange]-epsilon^2*splitRHS["TranspV*im_B*V"][Int(xrange/2)+1:xrange, 1:Int(xrange/2)]'))
=#


# check, if the sum over the rows adds up to 0
#=
dim_s = size(IBP)[1]
IBP_sumrow = zeros(dim_s)
for i in 1:dim_s
    for j in 1:dim_s 
        IBP_sumrow[i] += IBP[i, j] 
    end
end
display(IBP_sumrow)
=#
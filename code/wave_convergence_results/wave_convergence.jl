include("../functions.jl")
include("../analysis_tools_wave.jl")

using LinearAlgebra, Plots, PolynomialBases, DelimitedFiles, LaTeXStrings

for (ieps, epsilon) in enumerate([1.0 0.4 0.05])
    for basis in [GaussLegendre, LobattoLegendre]
        plot(xlabel = "steps", ylabel = "error", xscale = :log10, yscale = :log10)
        if basis == LobattoLegendre
            basissstr = "Lobatto"
        else
            basissstr = "Gauss"
        end	
        save_stepsizes = 0 # just possible to assign to next loop, because its constant in every iteration!!!!
        for deg in [0,1, 2, 3]
            CFLs = [1.0, 0.3, 0.02] # for the different epsilon, worked in practice
            CFL = CFLs[ieps] 
            stepsizes, errors = test_convergence_wave_background(LobattoLegendre, deg ; TMM = SSPRK10_4, exprange = [3,10], CFL = CFL*(1/(2*deg+1)), a=1, epsilon = epsilon)
            save_stepsizes = stepsizes # just possible to assign to next loop, because its constant in every iteration!!!!
            data = hcat(stepsizes, errors)
            writedlm(joinpath(@__DIR__,"./data/$(basissstr)_eps=$(epsilon)_order=$(deg+1)_data.txt"), data)
            plot!(stepsizes, errors, label = "order $(deg+1)")
            plot!(title = L"$\epsilon =  %$(epsilon)$, basis = %$(basissstr)")
        end
        savefig(joinpath(@__DIR__,"./plots/wave_$(basissstr)_eps=$(epsilon)_convergence.pdf"))

        for deg in [0,1, 2, 3]
            plot!(save_stepsizes, (save_stepsizes.+0.0).^(-(deg+1)), label = "reforder $(deg+1)")
        end
        savefig(joinpath(@__DIR__,"./plots/wave_$(basissstr)_eps=$(epsilon)_convergence(incl_reflines).pdf"))
    end
end




# read/plot scetch to split scripts maybe in future
#=
steps_and_error = readdlm(joinpath(@__DIR__,"./data/wave_convergence_results/Gauss_eps=0.05_order=2_data.txt"),'\t', Float64,'\n')
steps = steps_and_error[:,1]
error = steps_and_error[:,2]
plot(steps, error, xscale = :log10, yscale = :log10)
=#
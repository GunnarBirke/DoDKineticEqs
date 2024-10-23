using Plots
function vizualize(data, x_d, t_d; steplim = length(data[1,:]), u_exact = 0)  
    anim = @animate for i ∈ 1:steplim
        plot(legend = :outertopright, title = "T = $(t_d[i])")
        plot!(x_d, data[:,i], label = "u_approx", yrange = [-2.5, 2])
        if u_exact != 0
            plot!(x_d ,u_exact[:,i],label = "u_exact")
        end
    end
    return gif(anim, fps = 100)
end
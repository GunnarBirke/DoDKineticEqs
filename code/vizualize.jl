using Plots
function vizualize(data, x_d, t_d; steplim = min(length(data[1,:]), length(t_d)), u_exact = 0, clabel = 0) 
    if clabel == 0
        labelstring = ""
        else
        labelstring = clabel * ": "
    end
    ylimits = [minimum(minimum(data)), maximum(maximum(data))]
    anim = @animate for i ∈ 1:steplim
        plot(legend = :outertopright, title = "T = $(t_d[i])", ylims = ylimits)
        plot!(x_d, data[:,i], label = labelstring * "approx")
        if u_exact != 0
            plot!(x_d ,u_exact[:,i],label = labelstring * "exact", linestyle = :dashdotdot)
        end
    end
    return gif(anim, fps = 250)
end
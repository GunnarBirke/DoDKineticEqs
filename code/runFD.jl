include("Telegraph.equation.v1.0.jl")
include("vizualize.jl")


# solving
params = (D1=D1, ε=1.01);

tsteps = 50;
prob = ODEProblem(telegraph!, u0, tspan);
sol = solve(prob, Kvaerno3(), saveat=range(first(tspan), stop=last(tspan), length=tsteps));

# prepraring plot
t_vec=range(first(tspan), last(tspan), length = tsteps)
x_vec=range(x_min, x_max, length = N)
full_sol = sol[1:N, :] .+sol[N+1:2*N, :]*params.ε
plot(camera = (80, 30), xlabel = "t", ylabel = "x")

# plots:
plot(camera = (50, 35))
# plot f
plot!(t_vec, x_vec, full_sol,  st=:surface)
# plot rho
#plot!(t_vec, x_vec, sol[1:N, 1:tsteps],  st=:surface)
# plot g
#plot!(t_vec, x_vec, sol[N+1:2*N, 1:tsteps],  st=:surface)

vizualize(full_sol, x_vec, t_vec, clabel = "full")
#vizualize(sol[1:N, 1:tsteps], x_vec, t_vec, clabel = "rho")
#vizualize(sol[N+1:2*N, 1:tsteps], x_vec, t_vec, clabel = "g")
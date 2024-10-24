##############################
#####                    #####
##### Telegraph-equation #####
##### Version: 1.0       #####
#####                    #####
##############################

using SummationByPartsOperators
using LinearAlgebra
using OrdinaryDiffEq
using LaTeXStrings
using Plots: Plots, plot, plot!, savefig

########## initial data: square-shaped function ##########

# u0 = f(u,v) = (rho(x)=1,j=0), x ∈ [-0.2,0.2]; 0 elsewhere
# more details in: https://epubs.siam.org/doi/abs/10.1137/07069479X

function u1_func(x)
    if x >= -0.2 && x <= 0.2
        rho = 1.0
    else
        rho = 0.0
    end
end

function u2_func(x)
    j = 0.0
end

function telegraph!(du, u, p, t)
    u1 = u[1:N]
    u2 = u[N+1:2*N]
    
    du[1:N] = -params.D1 * u2                                                           # du[1:N] = du1
    du[N+1:2*N] = ((-1.0/(params.ε)^2) * (params.D1 * u1)) - (1.0/(params.ε)^2)*(u2)    # du[N+1:2*N] = du2
end

x_min   = -1.0;
x_max   = 1.0;
N       = 100;          # cell discretization, grid-points
tspan   = (0.0, 0.1);

D1 = periodic_derivative_operator(derivative_order=1, accuracy_order=4,
                                  xmin=x_min, xmax=x_max, N=N)

u0 = zeros(1:2*N);
u0[1:N] = u1_func.(grid(D1));
u0[N+1:2*N] = u2_func.(grid(D1));
params = (D1=D1, ε=1.0);


prob = ODEProblem(telegraph!, u0, tspan);
sol = solve(prob, Kvaerno3(), saveat=range(first(tspan), stop=last(tspan), length=5));
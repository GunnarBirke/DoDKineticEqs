include("Telegraph.equation.v1.0.jl")

D1 = periodic_derivative_operator(derivative_order=1, accuracy_order=1,
                                  xmin=x_min, xmax=x_max, N=N)

params = (D1=D1, ε=1.0);
println(D1)
function telegraph(u)
    du = zeros(2*N)
    u1 = u[1:N]
    u2 = u[N+1:2*N]
    
    du[1:N] = -params.D1 * u2                                                           # du[1:N] = du1
    du[N+1:2*N] = ((-1.0/(params.ε)^2) * (params.D1 * u1)) - (1.0/(params.ε)^2)*(u2)    # du[N+1:2*N] = du2
    return du
end


function get_telegraph_RHS(f,N)
    RHS_mat_tel = zeros(2*N,2*N)
    for i in 1:2*N
        unitvec = zeros(2*N)
        unitvec[i] = 1.0
        RHS_mat_tel[1:2*N, i] = telegraph(unitvec)
    end
    return RHS_mat_tel
end

RHS_mat_tel = get_telegraph_RHS(telegraph, 100)

RHS_mat_tel[101:110, 1:10]
RHS_mat_tel[1:10, 101:110]
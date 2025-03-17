function SSPRK3(setup, RHS_mat)
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    u_exact = get_exact_solution(setup)
    Nt = length(t_d)
    # Nx = length(x_d)*2, because wave/telegraph is a 2-dim system of equations
    Nx = length(x_d)*2
    u_sol = zeros(Nx, Nt)
    substeps = zeros(Nx, 3)
    f(u,t) = f_RHS(u,t, RHS_mat, setup)
    cRK = [0, 1, 1/2]
    for it in range(1,Nt, step = 1)
        if it == 1
            u_sol[:,it] = u0[1:Nx]
        else
            dt = (t_d[it]-t_d[it-1])
            substeps[:, 1] = f(u_sol[:,it - 1],t_d[it - 1] + cRK[1]*dt)
            substeps[:, 2] = f(u_sol[:,it - 1] + dt * (substeps[:, 1]), t_d[it - 1] + cRK[2]*dt)
            substeps[:, 3] = f(u_sol[:,it - 1] + dt * 1/4*(substeps[:, 1] + substeps[:, 2]), t_d[it - 1] + cRK[3]*dt)
            u_sol[:,it] = u_sol[:,it - 1] + dt * (1/6*substeps[:, 1] + 1/6*substeps[:, 2] + 2/3 * substeps[:, 3])
        end
    end
    return u_sol, u_exact
end

function SSPRK10_4(setup, RHS_mat)
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    u_exact = get_exact_solution(setup)
    Nt = length(t_d)
    Nx = length(x_d)*2
    u_sol = zeros(Nx, Nt)
    substeps = zeros(Nx, 10)
    f(u,t) = f_RHS(u,t, RHS_mat, setup)
    cRK = [0, 1/6, 1/3, 1/2, 2/3, 1/3, 1/2, 2/3, 5/6, 1]
    for it in range(1,Nt, step = 1)
        if it == 1
            u_sol[:,it] = u0[1:Nx]
        else
            dt = (t_d[it]-t_d[it-1])
            substeps[:, 1] = f(u_sol[:,it - 1], t_d[it - 1] + cRK[1]*dt)
            substeps[:, 2] = f(u_sol[:,it - 1] + dt * 1/6*(substeps[:, 1]), t_d[it - 1] + cRK[2]*dt)
            substeps[:, 3] = f(u_sol[:,it - 1] + dt * 1/6*(substeps[:, 1] + substeps[:, 2]), t_d[it - 1] + cRK[3]*dt)
            substeps[:, 4] = f(u_sol[:,it - 1] + dt * 1/6*(substeps[:, 1] + substeps[:, 2] + substeps[:, 3]), t_d[it - 1] + cRK[4]*dt)
            substeps[:, 5] = f(u_sol[:,it - 1] + dt * 1/6*(substeps[:, 1] + substeps[:, 2] + substeps[:, 3] + substeps[:, 4]), t_d[it - 1] + cRK[5]*dt)
            mid_sum15 = 1/15*(substeps[:, 1] + substeps[:, 2] + substeps[:, 3] + substeps[:, 4] + substeps[:, 5])
            substeps[:, 6] = f(u_sol[:,it - 1] + dt * mid_sum15, t_d[it - 1] + cRK[6]*dt)
            substeps[:, 7] = f(u_sol[:,it - 1] + dt * (mid_sum15 + 1/6*(substeps[:, 6])), t_d[it - 1] + cRK[7]*dt)
            substeps[:, 8] = f(u_sol[:,it - 1] + dt * (mid_sum15 + 1/6*(substeps[:, 6] + substeps[:, 7])), t_d[it - 1] + cRK[8]*dt)
            substeps[:, 9] = f(u_sol[:,it - 1] + dt * (mid_sum15 + 1/6*(substeps[:, 6] + substeps[:, 7] + substeps[:, 8])), t_d[it - 1] + cRK[9]*dt)
            substeps[:, 10] = f(u_sol[:,it - 1] + dt * (mid_sum15 + 1/6*(substeps[:, 6] + substeps[:, 7] + substeps[:, 8] + substeps[:, 9])), t_d[it - 1] + cRK[10]*dt)
            subsum = zeros(Nx)
            for isubsum in 1:10
                subsum += substeps[:, isubsum]
            end
            u_sol[:,it] = u_sol[:,it - 1] + dt * 1/10 * subsum
        end
    end
    return u_sol, u_exact
end
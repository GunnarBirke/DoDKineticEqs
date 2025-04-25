function ImExEuler(setup, RHS_mat)
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = problem["ex_RHS_mat"]
    im_RHS_mat = problem["im_RHS_mat"]
    # Instead of Nx, use the following notation (where dim = Nx for system, compdim = Nx for scalar equations)
    dim = size(ex_RHS_mat/2)[1]
    compdim = Int(dim/2)
    # Checking, if the splitting is possible
    im_RHS_mat[1:compdim, :] == zeros(compdim, dim) || throw(ArgumentError("pseudo-ImEx splitting not possible: implicit component has unexpected entrys "))
    # General Initialisation
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(dim, Nt)
    # Pseudo-ImEx specialisation
    u_sol_rho = zeros(compdim, Nt)
    u_sol_j = zeros(compdim, Nt)
    ex_RHS_j_to_rho = ex_RHS_mat[1:compdim, compdim + 1:dim]
    ex_RHS_j_to_j = ex_RHS_mat[compdim + 1:dim, compdim + 1:dim]
    im_RHS = im_RHS_mat[compdim + 1:dim, :]
    #println(cond((I-(t_d[2]-t_d[1])*im_RHS[:, compdim + 1:dim])))
    # as we restrict ourselves to periodic boundary conditions, calling f_RHS is not needed
    for it in range(1,Nt, step = 1)
        if it == 1
            u_sol_rho[:,it] = u0[1:compdim, it]
            u_sol_j[:,it] = u0[compdim + 1:dim, it]
            u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
        else
            dt = (t_d[it]-t_d[it-1])
            u_sol_rho[:, it] = u_sol_rho[:, it - 1] + dt*ex_RHS_j_to_rho*u_sol_j[:, it - 1]
            u_sol_j[:, it] = (I-dt*im_RHS[:, compdim + 1:dim])\(u_sol_j[:, it - 1] + dt*ex_RHS_j_to_j*u_sol_j[:, it - 1] + dt*im_RHS[:, 1:compdim]*u_sol_rho[:, it])
        end
        u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
    end

    return u_sol, u_exact
end

function ARS2(setup, RHS_mat)
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = problem["ex_RHS_mat"]
    im_RHS_mat = problem["im_RHS_mat"]
    # Instead of Nx, use the following notation (where dim = Nx for system, compdim = Nx for scalar equations)
    dim = size(ex_RHS_mat/2)[1]
    compdim = Int(dim/2)
    # Checking, if the splitting is possible
    im_RHS_mat[1:compdim, :] == zeros(compdim, dim) || throw(ArgumentError("pseudo-ImEx splitting not possible: implicit component has unexpected entrys "))
    # General Initialisation
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(dim, Nt)
    # Pseudo-ImEx specialisation
    u_sol_rho = zeros(compdim, Nt)
    u_sol_j = zeros(compdim, Nt)
    ex_RHS_j_to_rho = ex_RHS_mat[1:compdim, compdim + 1:dim]
    ex_RHS_j_to_j = ex_RHS_mat[compdim + 1:dim, compdim + 1:dim]
    im_RHS = im_RHS_mat[compdim + 1:dim, :]
    #println(cond((I-(t_d[2]-t_d[1])*im_RHS[:, compdim + 1:dim])))
    # as we restrict ourselves to periodic boundary conditions, calling f_RHS is not needed
    for it in range(1,Nt, step = 1)
        if it == 1
            u_sol_rho[:,it] = u0[1:compdim, it]
            u_sol_j[:,it] = u0[compdim + 1:dim, it]
            u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
        else
            dt = (t_d[it]-t_d[it-1])
            gamma = 1-1/sqrt(2)
            delta = 1-1/(2*gamma)
            # k_1
            k1_rho = u_sol_rho[:,it - 1]
            k1_j = u_sol_j[:,it - 1]
            # k_2
            k2_rho = u_sol_rho[:,it - 1] + dt*gamma*ex_RHS_j_to_rho*k1_j
            k2_j = (I-dt*gamma*im_RHS[:, compdim + 1:dim])\(u_sol_j[:, it - 1] + dt*gamma*ex_RHS_j_to_j*k1_j + dt*gamma*im_RHS[:, 1:compdim]*k2_rho)
            # k_3
            k3_rho =  u_sol_rho[:,it - 1] + dt*(delta*ex_RHS_j_to_rho*k1_j + (1-delta)*ex_RHS_j_to_rho*k2_j)
            k3_j = (I-dt*gamma*im_RHS[:, compdim + 1:dim])\(u_sol_j[:, it - 1] + dt*(delta*ex_RHS_j_to_j*k1_j + (1-delta)*ex_RHS_j_to_j*k2_j + (1-gamma)*im_RHS*vcat(k2_rho, k2_j) + gamma*im_RHS[:, 1:compdim]*k3_rho))
            # aufgrund von gobally stiffly accurate:
            u_sol_rho[:, it] = k3_rho
            u_sol_j[:, it] = k3_j
            # ganze pseudo-IMEX-Prozedur kann man sicher auch noch etwas allgemeiner aufziehen, ohne es wirklich komplizieter zu machen: Hier haben wir die drei Teile 
            # ex_RHS_j_to_rho, x_RHS_j_to_j, im_RHS[:, 1:compdim], im_RHS[:, compdim + 1:dim]. Allgemein kann man diese 4 Teile sicher besser identifizieren und herausarbeiten,
            # sodass dies alles übersichlicher und allgemeiner wird. Hier jetzt aber erstmal Stumpf für unser System.
            # beachte auch, dass der explizite Teil nur so speziell auftritt, da er an vielen Stellen 0 ist. Ansonste würden mehr explizite Einträge sowohl in rho,
            # als auch j vorkommen (gleiches Argument gilt für impliziten Teil, wobei dieser für die pseudo-IMEX betrachtung nicht großartig anders aussehen darf)
        end
        u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
    end

    return u_sol, u_exact
end

function ARS3(setup, RHS_mat)
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = problem["ex_RHS_mat"]
    im_RHS_mat = problem["im_RHS_mat"]
    # Instead of Nx, use the following notation (where dim = Nx for system, compdim = Nx for scalar equations)
    dim = size(ex_RHS_mat/2)[1]
    compdim = Int(dim/2)
    # Checking, if the splitting is possible
    im_RHS_mat[1:compdim, :] == zeros(compdim, dim) || throw(ArgumentError("pseudo-ImEx splitting not possible: implicit component has unexpected entrys "))
    # General Initialisation
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(dim, Nt)
    # Pseudo-ImEx specialisation
    u_sol_rho = zeros(compdim, Nt)
    u_sol_j = zeros(compdim, Nt)
    ex_RHS_j_to_rho = ex_RHS_mat[1:compdim, compdim + 1:dim]
    ex_RHS_j_to_j = ex_RHS_mat[compdim + 1:dim, compdim + 1:dim]
    im_RHS = im_RHS_mat[compdim + 1:dim, :]
    #println(cond((I-(t_d[2]-t_d[1])*im_RHS[:, compdim + 1:dim])))
    # as we restrict ourselves to periodic boundary conditions, calling f_RHS is not needed
    for it in range(1,Nt, step = 1)
        if it == 1
            u_sol_rho[:,it] = u0[1:compdim, it]
            u_sol_j[:,it] = u0[compdim + 1:dim, it]
            u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
        else
            dt = (t_d[it]-t_d[it-1])
            # k_1
            k1_rho = u_sol_rho[:,it - 1]
            k1_j = u_sol_j[:,it - 1]
            # k_2
            k2_rho = u_sol_rho[:,it - 1] + dt*1/2*ex_RHS_j_to_rho*k1_j
            k2_j = (I-dt*1/2*im_RHS[:, compdim + 1:dim])\(u_sol_j[:, it - 1] + dt*1/2*ex_RHS_j_to_j*k1_j + dt*1/2*im_RHS[:, 1:compdim]*k2_rho)
            # k_3
            k3_rho = u_sol_rho[:,it - 1] + dt*(11/18*ex_RHS_j_to_rho*k1_j + 1/18*ex_RHS_j_to_rho*k2_j)
            k3_j = (I-dt*1/2*im_RHS[:, compdim + 1:dim])\(u_sol_j[:, it - 1] + dt*(11/18*ex_RHS_j_to_j*k1_j + 1/18*ex_RHS_j_to_j*k2_j + 1/6*im_RHS*vcat(k2_rho, k2_j) + 1/2*im_RHS[:, 1:compdim]*k3_rho))
            # k_4
            k4_rho = u_sol_rho[:,it - 1] + dt*(ex_RHS_j_to_rho*(5/6*k1_j + -5/6*k2_j + 1/2*k3_j))
            k4_j = (I-dt*1/2*im_RHS[:, compdim + 1:dim])\(u_sol_j[:, it - 1] + dt*(ex_RHS_j_to_j*(5/6*k1_j + -5/6*k2_j + 1/2*k3_j)
                     + im_RHS*(-1/2*vcat(k2_rho, k2_j) + 1/2*vcat(k3_rho, k3_j)) + 1/2*im_RHS[:, 1:compdim]*k4_rho))
            # k_5
            k5_rho = u_sol_rho[:,it - 1] + dt*(ex_RHS_j_to_rho*(1/4*k1_j + 7/4*k2_j + 3/4*k3_j + -7/4*k4_j))
            k5_j = (I-dt*1/2*im_RHS[:, compdim + 1:dim])\(u_sol_j[:, it - 1] + dt*(ex_RHS_j_to_j*(1/4*k1_j + 7/4*k2_j + 3/4*k3_j + -7/4*k4_j) 
                     +  im_RHS*(3/2*vcat(k2_rho, k2_j) + -3/2*vcat(k3_rho, k3_j) + 1/2*vcat(k4_rho, k4_j)) + 1/2*im_RHS[:, 1:compdim]*k5_rho))
            # aufgrund von gobally stiffly accurate:
            u_sol_rho[:, it] = k5_rho
            u_sol_j[:, it] = k5_j
        end
        u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
    end

    return u_sol, u_exact
end

function SSPRK3(setup, RHS_mat)
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
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
    #u_exact = 0
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
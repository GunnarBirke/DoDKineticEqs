function ImExEuler(setup, RHS_mat)
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = setup["ex_RHS_mat"]
    im_RHS_mat = setup["im_RHS_mat"]
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
    eq_type = setup["eq_type"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(dim, Nt)
    if eq_type == "telegraph"
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
    elseif eq_type == "heat"
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol[:, it] = u0[:, it]
            else
                dt = (t_d[it]-t_d[it-1])
                u_sol[:, it] = u_sol[:, it - 1] + dt*RHS_mat*u_sol[:, it - 1]
            end
        end
    end
    return u_sol, u_exact
end

# treating also the A(also denoted by b_h) implicit for the ImEx Euler method
function ImExEuler2(setup, RHS_mat)
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = setup["ex_RHS_mat"]
    im_RHS_mat = setup["im_RHS_mat"]
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
    eq_type = setup["eq_type"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(dim, Nt)
    if eq_type == "telegraph"
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
                u_sol_j[:, it] = (I-dt*im_RHS[:, compdim + 1:dim]-dt*ex_RHS_j_to_j)\(u_sol_j[:, it - 1] + dt*im_RHS[:, 1:compdim]*u_sol_rho[:, it])
            end
            u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
        end
    elseif eq_type == "heat"
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol[:, it] = u0[:, it]
            else
                dt = (t_d[it]-t_d[it-1])
                u_sol[:, it] = u_sol[:, it - 1] + dt*RHS_mat*u_sol[:, it - 1]
            end
        end
    end
    return u_sol, u_exact
end

function ARS2(setup, RHS_mat)
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = setup["ex_RHS_mat"]
    im_RHS_mat = setup["im_RHS_mat"]
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
    eq_type = setup["eq_type"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(dim, Nt)
    if eq_type == "telegraph"
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
    elseif eq_type == "heat"
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol[:, it] = u0[:, it]
            else
                dt = (t_d[it]-t_d[it-1])
                gamma = 1-1/sqrt(2)
                delta = 1-1/(2*gamma)
                k1 = u_sol[:, it-1]
                k2 = u_sol[:,it - 1] + dt*gamma*RHS_mat*k1
                k3 = u_sol[:,it - 1] + dt*(delta*RHS_mat*k1 + (1-delta)*RHS_mat*k2)
                u_sol[:, it] = k3
            end
        end
    end
    return u_sol, u_exact
end

function ARS3(setup, RHS_mat)
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = setup["ex_RHS_mat"]
    im_RHS_mat = setup["im_RHS_mat"]
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
    eq_type = setup["eq_type"]
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
    if eq_type == "telegraph"
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
    elseif eq_type == "heat"
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol[:, it] = u0[:, it]
            else
                dt = (t_d[it]-t_d[it-1])
                k1 = u_sol[:, it-1]
                k2 = u_sol[:,it - 1] + dt*1/2*RHS_mat*k1
                k3 = u_sol[:,it - 1] + dt*(11/18*RHSmat*k1 + 1/18*RHS_mat*k2)
                k4 = u_sol[:,it - 1] + dt*(RHSmat*(5/6*k1 + -5/6*k2 + 1/2*k3))
                k5 = u_sol[:,it - 1] + dt*(RHSmat*(1/4*k1 + 7/4*k2 + 3/4*k3 + -7/4*k4))
                u_sol[:, it] = k5
            end
        end
    end
    return u_sol, u_exact
end

function ImEx(setup, RHS_mat, Tableau; only_explicit = false)
    # matrices
    A_ex = Tableau["A_ex"]
    A_im = Tableau["A_im"]
    b_ex = Tableau["b_ex"]
    b_im = Tableau["b_im"]
    s = length(b_im)
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = setup["ex_RHS_mat"]
    im_RHS_mat = setup["im_RHS_mat"]
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
    eq_type = setup["eq_type"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(dim, Nt)
    if (eq_type in ["telegraph", "telegraph_symm", "wave", "wave_symm"]) && only_explicit == false
        # Pseudo-ImEx specialisation
        u_sol_rho = zeros(compdim, Nt)
        u_sol_j = zeros(compdim, Nt)
        ex_RHS_j_to_rho = ex_RHS_mat[1:compdim, compdim + 1:dim]
        ex_RHS_j_to_j = ex_RHS_mat[compdim + 1:dim, compdim + 1:dim]
        im_RHS = im_RHS_mat[compdim + 1:dim, :]
       
        # as we restrict ourselves to periodic boundary conditions, calling f_RHS is not needed
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol_rho[:,it] = u0[1:compdim, it]
                u_sol_j[:,it] = u0[compdim + 1:dim, it]
                u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
            else
                dt = (t_d[it]-t_d[it-1])
                k_rho = zeros(compdim, s)
                k_j = zeros(compdim, s)
                for i in 1:s
                    k_rho[:, i] = u_sol_rho[:,it - 1]
                    k_j[:, i] = u_sol_j[:,it - 1]
                    for l = 1:i-1
                        k_rho[:, i] += dt*A_ex[i, l]*ex_RHS_j_to_rho*k_j[:, l]
                        k_j[:, i] += dt*A_ex[i, l]*ex_RHS_j_to_j*k_j[:, l] + dt*A_im[i, l]*im_RHS*vcat(k_rho[:, l], k_j[:, l])
                    end
                    k_j[:, i] += dt*A_im[i,i]*im_RHS[:, 1:compdim]*k_rho[:, i]
                    k_j[:, i] = (I-dt*A_im[i,i]*im_RHS[:, compdim + 1:dim])\k_j[:, i]
                end
                # TODO
                # Dazu im besten Fall noch Checks einabuen, dass A_ex wirklich explizit und A_im wirklich lower triangular ist
                u_sol_rho[:, it] = k_rho[:, s]
                u_sol_j[:, it] = k_j[:, s]
                
                u_sol_rho[:, it] = u_sol_rho[:, it - 1]
                u_sol_j[:, it] = u_sol_j[:, it - 1]
                for i in 1:s
                    u_sol_rho[:, it] += dt*b_ex[i]*ex_RHS_j_to_rho*k_j[:, i]
                    u_sol_j[:, it] += dt*b_ex[i]*ex_RHS_j_to_j*k_j[:, i] + dt*b_im[i]*im_RHS*vcat(k_rho[:, i], k_j[:, i])
                end
                
            end
            u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
        end
    elseif (eq_type in ["heat", "transport"]) || only_explicit == true
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol[:, it] = u0[:, it]
            else
                dt = (t_d[it]-t_d[it-1])
                k = zeros(dim, s)
                for i in 1:s
                    k[:, i] = u_sol[:,it - 1]
                    for l = 1:i-1
                        k[:, i] += dt*A_ex[i, l]*RHS_mat*k[:, l]
                    end
                end
                u_sol[:, it] = u_sol[:,it - 1]
                for i in 1:s
                    u_sol[:,it] += dt*b_ex[i]*RHS_mat*k[:, i]
                end
            end
        end
    end
    return u_sol, u_exact
end


function get_RK_tableau(tableau_string)
    Tableau = Dict()
    if tableau_string == "ImExEuler" #FSAL, SA, GSA
        Tableau["A_ex"] = [0 0; 1 0]
        Tableau["A_im"] = [0 0; 0 1]
        Tableau["b_ex"] = [1 0]
        Tableau["b_im"] = [0 1]
    elseif tableau_string == "ARS2" #FSAL, SA, GSA
        gamma = 1-1/sqrt(2)
        delta = 1-1/(2*gamma)
        Tableau["A_ex"] = [0 0 0; gamma 0 0; delta 1-delta 0]
        Tableau["A_im"] = [0 0 0; 0 gamma 0; 0 1-gamma gamma]
        Tableau["b_ex"] = [delta 1-delta 0]
        Tableau["b_im"] = [0 1-gamma gamma]
    elseif tableau_string == "ARS3" #FSAL, SA, GSA
        Tableau["A_ex"] = [0 0 0 0 0;
        1/2 0 0 0 0;
        11/18 1/18 0 0 0;
        5/6 -5/6 1/2 0 0;
        1/4 7/4 3/4 -7/4 0]
        Tableau["A_im"] = [0 0 0 0 0;
                    0 1/2 0 0 0;
                    0 1/6 1/2 0 0;
                    0 -1/2 1/2 1/2 0;
                    0 3/2 -3/2 1/2 1/2]
        Tableau["b_ex"] = [1/4 7/4 3/4 -7/4 0]
        Tableau["b_im"] = [0, 3/2, -3/2, 1/2, 1/2]
    elseif tableau_string == "AGSA342" #FSAL, SA, GSA (order 2, Type I ImEx)
        Tableau["A_ex"] = [ 0 0 0 0;
                    -139833537/38613965 0 0 0;
                    85870407/49798258 -121251843/1756367063 0 0;
                    1/6 1/6 2*1/3 0]
        Tableau["A_im"] = [ 168999711/74248304 0 0 0;
                    44004295/24775207 202439144/118586105 0 0;
                    -6418119/169001713 -748951821/1043823139 12015439/183058594 0;
                    -370145222/355758315 1/3 0 202439144/118586105]
        Tableau["b_ex"] = [1/6 1/6 2/3 0]
        Tableau["b_im"] = [-370145222/355758315 1/3 0 202439144*1/118586105]
    elseif tableau_string == "IGSA2" #FSAL, SA, GSA (order 2, Type I ImEx)
        Tableau["A_ex"] = [ 0 0 0 0;
                    1/3 0 0 0;
                    7/24 3/8 0 0;
                    1/2 -1/2 1 0]
        Tableau["A_im"] = [ 1/4 0 0 0;
                    0 1/4 0 0;
                    1/16 3/16 1/4 0;
                    1/4 1/4 1/4 1/4]
        Tableau["b_ex"] = [1/2 -1/2 1 0]
        Tableau["b_im"] = [1/4 1/4 1/4 1/4]
    elseif tableau_string == "SSP2ImEx332" # not FSAL, SA, not GSA (order 2, Type I ImEx) 
        Tableau["A_im"] = [1/4 0 0;
                    0 1/4 0;
                    1/3 1/3 1/3]
        Tableau["b_im"] = [1/3 1/3 1/3]
        Tableau["A_ex"] = [0 0 0;
                    1/2 0 0;
                    1/2 1/2 0]
        Tableau["b_ex"] = [1/3 1/3 1/3]
    elseif tableau_string == "GSA1" # FSAL, SA, GSA (order 1, Type I ImEx) 
        gamma = 0.5 # free choice for gamma > 0
        a_ex = 0.5 # free choice for gamma > 0
        w1_ex = 0.5 # free choice for gamma > 0
        Tableau["A_im"] = [gamma 0 0;
                    0 gamma 0;
                    1-gamma 0 gamma]
        Tableau["b_im"] = [1-gamma 0 gamma]
        Tableau["A_ex"] = [0 0 0;
                    a_ex 0 0;
                    w1_ex 1-w1_ex 0]
        Tableau["b_ex"] = [w1_ex 1-w1_ex 0]
    elseif tableau_string == "SSP3ImEx343" # not FSAL, not SA, not GSA (order 3, Type I ImEx)
        α = 0.24169426078821
        β = 0.06042356519705
        η = 0.12915286960590
        Tableau["A_im"] = [ α 0 0 0;
                   -α α 0 0;
                    0 1-α α 0;
                   β η 1/2−β−η−α α]
        Tableau["b_im"] = [0 1/6 1/6 2/3]
    
        Tableau["A_ex"] = [0 0 0 0;
                      0 0 0 0;
                      0 1 0 0;
                      0 1/4 1/4 0]
        Tableau["b_ex"] = [0 1/6 1/6 2/3]
    elseif tableau_string == "BPR343" # No Type 1, but FSAL, SA
        Tableau["A_im"] = [0 0 0 0 0;
                    1/2 1/2 0 0 0;
                    5/18 -1/9 1/2 0 0;
                    1/2 0 0 1/2 0;
                    1/4 0 3*1/4 -1/2 1/2]
        Tableau["b_im"] = [1/4, 0, 3/4, -1/2, 1/2]
        Tableau["A_ex"] = [0 0 0 0 0;
                    1 0 0 0 0;
                    4/9 2/9 0 0 0;
                    1/4 0 3/4 0 0;
                    1/4 0 3/4 0 0]
        Tableau["b_ex"] = [1/4, 0, 3/4, 0, 0]
    elseif tableau_string == "ImExEulernotGSA"
        Tableau["A_ex"] = [0]
        Tableau["A_im"] = [1]
        Tableau["b_ex"] = [1]
        Tableau["b_im"] = [1]
    end
    return Tableau
end


function Heun(setup, RHS_mat)
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    eq_type = setup["eq_type"]
    u_exact = get_exact_solution(setup)
    Nt = length(t_d)
    Nx = length(u0)
    u_sol = zeros(Nx, Nt)
    substeps = zeros(Nx, 2)
    f(u,t) = f_RHS(u,t, RHS_mat, setup)
    cRK = [0, 1]
    for it in range(1,Nt, step = 1)
        if it == 1
            u_sol[:,it] = u0[1:Nx]
        else
            dt = (t_d[it]-t_d[it-1])
            substeps[:, 1] = f(u_sol[:,it - 1],t_d[it - 1] + cRK[1]*dt)
            substeps[:, 2] = f(u_sol[:,it - 1] + dt * (substeps[:, 1]), t_d[it - 1] + cRK[2]*dt)
            u_sol[:,it] = u_sol[:,it - 1] + dt * (1/2*substeps[:, 1] + 1/2*substeps[:, 2])
        end
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
    Nx = length(u0)
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
    Nx = length(u0)
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
function setup_problem_eq(init_cond, xmin, xmax, cellnumber; Tmax = 1, a = 1, CFL = 0.9, bcs = "periodic", epsilon = 1.0)
    vertices = Array(range(xmin, xmax, length = cellnumber + 1))
    h = vertices[2]-vertices[1]
    t_d = range(0,Tmax, length = Integer(round(abs(Tmax/(CFL*h/a)))))
    # Setting initial condition as function
    random_ic = zeros(3*cellnumber)
    if init_cond == "random"
        random_ic = randn(3*cellnumber)
    end
    function u0_eval(x)
        if init_cond == "sin"        
            return sin(2*pi*x)
        elseif init_cond == "cos"
            return cos(4*pi*x)
        elseif init_cond == "squared"
            return (x-1/2)^2
        elseif init_cond == "const"
            return 1.0
        elseif init_cond == "linear"
            return x
        elseif init_cond == "hat"
            return hat(x)
        elseif init_cond == "discont"
            if x<1/2
                return -1.0
            else
                return 1.0
            end
        elseif init_cond == "plateau"
            if x >= 0.45 && x <= 0.55
                rho = 1.0
            else
                rho = 0.0
            end
        elseif init_cond == "telsin"
            r = -2/(1+sqrt(1-4*epsilon^2))
            return 1/r*sin(x)
        elseif init_cond == "random"
            return random_ic[ceil(Int, abs(x/(xmax-xmin))*cellnumber)+1]
        end
    end
    
    problem = Dict("a" => a, "u0_eval" => u0_eval, "v" => vertices, "Tmax" => Tmax, "h" => h, "CFL" => CFL,
                    "cellnumber" => cellnumber, "t_d" => t_d, "init_cond" => init_cond, "xmin" => xmin, "xmax" => xmax, "bcs" => bcs, "epsilon" => epsilon)
    return problem
end

function hat(x)
    if x > 0.4 && x < 0.6
        return 0.5-abs(x-0.5)
    else
        return 0.4
    end
end

function to_periodic(x)
    if x > 1 
        return x%1
    elseif x<0
        while x<0
            x = x+1
        end
        return x
    else 
        return x
    end
end

function get_exact_solution(setup)
    eq_type = setup["eq_type"]
    u0_eval = setup["u0_eval"]
    x_d = setup["x_d"]
    a = setup["a"]
    t_d = setup["t_d"]
    N = setup["cellnumber"]
    deg = setup["deg"]
    epsilon = setup["epsilon"]
    if setup["bcs"] == "periodic"
        #return u_exact = [u0_eval(to_periodic(x-(a*(t)))) for x in x_d, t in t_d]
        if eq_type == "wave"
            u_wave_exact = zeros((deg+1)*2*N, length(t_d));
            u_wave_exact[1:(deg+1)*N, :] = [1/2*(sin(2*pi*(x-1/epsilon*t))+epsilon*cos(2*pi*(x-1/epsilon*t))+sin(2*pi*(x+1/epsilon*t))-epsilon*cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in t_d];
            u_wave_exact[(deg+1)*N+1:(deg+1)*2*N, :] = [1/2*(1/epsilon*sin(2*pi*(x-1/epsilon*t))+cos(2*pi*(x-1/epsilon*t))-1/epsilon*sin(2*pi*(x+1/epsilon*t))+cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in t_d];
            return u_wave_exact
        elseif eq_type == "telegraph"
            println("not implemented yet")
        else 
            println("wrong eq type")
        end
    end
end

function f_RHS(u, t, RHS_mat, setup)
    u0_eval = setup["u0_eval"]
    x_d = setup["x_d"]
    a = setup["a"]
    basis = setup["basis"]
    deg = setup["deg"]
    v = setup["v"]
    # Nx = length(x_d)*2, because wave/telegraph is a 2-dim system of equations
    Nx = length(x_d)*2
    g_L = zeros(Nx)
    basis1 = basis
    for j in 1:deg+1
        basis1.nodes[j] = x_d[j] 
    end
    if setup["bcs"] == "dirichlet"
        Vs = interpolation_matrix(setup["xmin"], basis1)
        g_L[1:deg+1] = inv(diagm(basis.weights))*Vs'*u0_eval(-a*(t)) *2/(v[2]-v[1])
    end
    #(letzter Faktor ist invscale(1))
    #display(t)
    #display(u0_eval(-a*(t)))
    #display(g_L[1:6])
    return RHS_mat * u + a* g_L
end

function include_cut_cell(setup, cutcell_size, cutcell_pos)
    vert = setup["v"]
    create_cutcell = false
    if cutcell_size <=0.5
        if vert[cutcell_pos+1]-vert[cutcell_pos] > 1/2* setup["h"] && vert[cutcell_pos]-vert[cutcell_pos-1] > 1/2* setup["h"]
            create_cutcell = true
        end
    end
    if cutcell_size >0.5
        if vert[cutcell_pos+2]-vert[cutcell_pos + 1] > 1/2* setup["h"] && vert[cutcell_pos + 1]-vert[cutcell_pos] > 1/2* setup["h"]
            create_cutcell = true
        end
    end
    if create_cutcell
        vert_new = zeros(length(vert)+1)
        for i in range(1, length(vert)+1)
            if i < cutcell_pos
                vert_new[i] = vert[i] 
            elseif i == cutcell_pos
                vert_new[i] = vert[i] 
                vert_new[i+1] = ((1-cutcell_size)*vert[cutcell_pos] + (cutcell_size)*vert[cutcell_pos+1])
            elseif i > cutcell_pos + 1
                vert_new[i] = vert[i-1]
            end
        end
        setup["v"] = vert_new
        setup["cellnumber"] += 1
    else
        println("No new cut-cell created due to its neighbour being already a cut-cell")
    end
    return setup
end


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

####################################################################################################################################################
########################################################                          ##################################################################
########################################################    Telegraph equation    ##################################################################
########################################################                          ##################################################################
####################################################################################################################################################
function mymod(x, y)
    if x % y == 0
        return y
    else 
        return x % y
    end
end

function DGsemidiscretization_DoD_telegraph(setup, deg, nodes, num_flux; T = Float64, do_stabilize = false, c = 1.0, fix_eta = false, eq_type = "telegraph", fluxtype = "full")
    cellnumber = setup["cellnumber"]
    v = setup["v"]
    a = setup["a"]
    u0_eval = setup["u0_eval"]
    t_d = setup["t_d"]
    h = setup["h"]
    setup["nodes"] = nodes
    setup["deg"] = deg
    setup["num_flux"] = num_flux

    setup["eq_type"] = eq_type

    ### Telegraph parameters
    epsilon = setup["epsilon"]
    # fixed
    a = 1
    setup["bcs"] = "periodic"

    splitRHS = Dict()

    # basis defines the method, basis1 will be manipulated to serve interpolation matrices for other nodes
    basis = nodes(deg, T)
    basis1 = nodes(deg, T)

    # Some Parameters
    xrange = 2*cellnumber*(deg+1)
    function invscale(i, j)
        if j < xrange/2
            a*2/(v[i+1]-v[i])
        else
            2/(v[i+1]-v[i])
        end
    end
    shift(i) = (v[i+1]+v[i])/2
    scale(i) = (v[i+1]-v[i])/2
    
    # x_d, RHS calculation
    x_d = zeros(T, Int(xrange/2))
    M = zeros(T, xrange, xrange)
    M_global = zeros(T, xrange, xrange) # TODO
    DM = zeros(T, xrange, xrange)
    sM = zeros(T, xrange, xrange)
    DM_J1 = zeros(T, xrange, xrange)
    B = zeros(T, cellnumber*4, cellnumber*4)
    B_J0 = zeros(T, cellnumber*4, cellnumber*4)
    B_J0E1 = zeros(T, cellnumber*4, cellnumber*4)
    B_J0E2 = zeros(T, cellnumber*4, cellnumber*4)
    V = zeros(T, cellnumber*4, xrange)
    V_J0E1 = zeros(T, cellnumber*4, xrange)
    V_J0E2 = zeros(T, cellnumber*4, xrange)
    TranspV = zeros(T, xrange, cellnumber*4)
    Sc = zeros(T, xrange, xrange)
    
    # Creating Transformation Matrix from reference element to real elements
    for i in 1:2*cellnumber
        Sc[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] = diagm(ones(T, deg+1)*invscale(mymod(i, cellnumber),i))
    end
    # Constructing the interpolation matrices to the vertices of a cell at reference element
    left_interpolation = interpolation_matrix(-1 , basis1)
    right_interpolation = interpolation_matrix(1, basis1)
    for i in 1:cellnumber
        for j in 1:deg+1
            # Constructing the nodes
            x_d[(deg+1)*(i-1) + j] = basis.nodes[j]*scale(i) + shift(i)
            # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
            basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
        end
    end
    for i in 1:2*cellnumber
        # Constructing the Mass- and Differentiation Matrix
        M[(deg+1)*(i-1)+1:(deg+1)*i, (deg+1)*(i-1)+1:(deg+1)*i] = diagm(basis.weights)
        M_global[(deg+1)*(i-1)+1:(deg+1)*i, (deg+1)*(i-1)+1:(deg+1)*i] = diagm(basis.weights) * scale(mymod(i, cellnumber))
        # differentiate between rho and j, distribute correctly
        if i <= cellnumber
            DM[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(cellnumber+(i-1)) + 1:(deg+1)*(cellnumber+i)] = basis.D' * diagm(basis.weights)
        else
            DM[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*((i-1)-cellnumber) + 1:(deg+1)*(i-cellnumber)] = basis.D' * diagm(basis.weights)/(epsilon^2)
            sM[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] = diagm(basis.weights)/(epsilon^2)
        end

        # Constructing the Interpolation to the vertices for the flux terms with Gauss-Legendre nodes
        #V[2*i-1, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i], basis1)
        #V[2*i, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+1], basis1)
        V[2*i-1, (deg+1)*(i-1) + 1:(deg+1)*i] = left_interpolation
        V[2*i, (deg+1)*(i-1) + 1:(deg+1)*i] = right_interpolation
        
        # Extrapolation to the Inflow of Cell + 1(V2) and the Outflow of Cut-Cell(V3) from Cell - 1 (from the perspective of Cell -1).
        if do_stabilize
            if i < cellnumber
                if 2*(v[i+2]-v[i+1])/h<1
                    if fix_eta == false
                        eta = 1-minimum([1/(2*deg+1)*(v[i+1]-v[i])/((t_d[2]-t_d[1])*a), 1.0])
                    else
                        eta = 1-minimum([(v[i+1]-v[i])/h*1/c, 1.0])
                    end
                    V_J0E1[2*i+2, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+2], basis1) 
                    V_J0E2[2*i+3, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+2], basis1)
                end
            end
        end
            
        TranspV[(deg+1)*(i-1)+1:(deg+1)*i, 2*i-1:2*i] = V[2*i-1:2*i, (deg+1)*(i-1) + 1:(deg+1)*i]'
    end

    # Constructing numerical fluxes(with periodic BC)
    B = flux_terms_telegraph!(B, setup, deg, fluxtype = fluxtype)
    # DoD-Stabilization
    if do_stabilize
        for i in 1:cellnumber
            if 2*(v[i+1]-v[i])/h<1
                if fix_eta == false
                    eta = 1-minimum([1/(2*deg+1)*(v[i+1]-v[i])/((t_d[2]-t_d[1])*a), 1.0])
                else
                    eta = 1-minimum([(v[i+1]-v[i])/h*1/c, 1.0])
                end
                # Reduced Outflow  from Cut-Cell
                B_J0[2*i, 2*i] += -1 * eta
                # Reduced Inflow into Cell + 1 from Cut-Cell
                B_J0[2*(i) + 1, 2*(i)] += 1 * eta
                # Extended Outflowflow Cut-Cell from Cell -1 (corresponds to Reduced Inflow into Cut-Cell)
                B_J0E1[2*i , 2*i] += 1 * eta
                # Extended Inflow into Cell + 1 from Cell - 1
                B_J0E2[2*(i) + 1, 2*(i) + 1]+= -1 * eta
                
                # Volume Term
                if deg>0
                    DM_J1[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] +=
                        - basis.D' * diagm(basis.weights) * eta
                    
                    basis2=nodes(deg)
                    for j in 1:deg+1
                        # Inserting these nodes into the structure basis2, to adapt the field basis1.interp_matrix to the used nodes.
                        basis2.nodes[j] = x_d[(deg+1)*(i-2) + j] 
                    end
                    Vol_interpol = interpolation_matrix(x_d[(deg+1)*(i-1) + 1:(deg+1)*(i)] , basis2)
                    DM_J1[(deg+1)*(i-1) + 1:(deg+1)*(i), (deg+1)*(i-2) + 1:(deg+1)*(i-1)] += 
                        + basis.D' * diagm(basis.weights) * Vol_interpol * eta
                    # Extra-Stabilityterm with the extended test-function (L_E_in(u_h)-u_h)*L_E_in(v_h)
                    DM_J1[(deg+1)*(i-2) + 1:(deg+1)*(i-1), (deg+1)*(i-1) + 1:(deg+1)*i] +=
                        Vol_interpol' *  basis.D' * diagm(basis.weights) * eta
                    DM_J1[(deg+1)*(i-2) + 1:(deg+1)*(i-1), (deg+1)*(i-2) + 1:(deg+1)*(i-1)] += 
                        - Vol_interpol' * basis.D' * diagm(basis.weights) * Vol_interpol * eta
                end
            end
        end
    end
    
    # Constructing RHS Matrix
    invM = inv(M)
    if eq_type == "telegraph"
        RHS_mat = Sc*invM * (-TranspV* (B*V + B_J0*V + B_J0E1*V_J0E1 + B_J0E2*V_J0E2) + DM + DM_J1 - sM)
    elseif eq_type == "wave"
        #Wave equation (choosing epsilon =1/c)
        RHS_mat = Sc*invM * (-TranspV* (B*V + B_J0*V + B_J0E1*V_J0E1 + B_J0E2*V_J0E2) + DM + DM_J1)
    end


    # Store parts in dictionary
    splitRHS["invM"] = invM
    splitRHS["TranspV"] = TranspV
    splitRHS["B"] = B
    splitRHS["V"] = V
    splitRHS["B_J0"] = B_J0
    splitRHS["B_J0E1"] = B_J0E1
    splitRHS["B_J0E2"] = B_J0E2
    splitRHS["V_J0E1"] = V_J0E1
    splitRHS["V_J0E2"] = V_J0E2
    splitRHS["DM"] = DM
    splitRHS["DM_J1"] = DM_J1
    splitRHS["TransD"] = basis.D'
    splitRHS["M_used"] = M
    splitRHS["M_global"] = M_global
    splitRHS["Sc"] = Sc

    splitRHS["Flux_terms"] = Sc*invM * (-TranspV* (B*V + B_J0*V + B_J0E1*V_J0E1 + B_J0E2*V_J0E2))
    splitRHS["Flux_B"] = Sc*invM* -TranspV * B*V
    splitRHS["Flux_B_J0"] = Sc*invM* -TranspV * B_J0*V
    splitRHS["Flux_B_J0E1"] = Sc*invM* -TranspV * B_J0E1*V_J0E1
    splitRHS["Flux_B_J0E2"] = Sc*invM* -TranspV * B_J0E2*V_J0E2
    splitRHS["Volume_terms"] = Sc*invM * (DM + DM_J1)
    # Constructing the initial conition (special case for j coordinate)
    u0_rho = [u0_eval(x) for x in x_d]
    if eq_type == "telegraph"
        u0_j = cos.(x_d)
    elseif eq_type == "wave"
        u0_j = cos.(2*pi*x_d)
    end

    u0 = vcat(u0_rho, u0_j);
    setup["u0"] = u0
    setup["M"] = M_global
    setup["x_d"] = x_d
    setup["basis"] = basis
    #display((-TranspV* (B*V + B2*V2))[1:9, 1:9])
    #display(DM_J1[1:10, 1:10])
    
    return RHS_mat, setup, splitRHS
end


function flux_terms_telegraph!(B, setup, deg; fluxtype = "altlr")
    cellnumber = setup["cellnumber"]
    epsilon = setup["epsilon"]
    eps2 = 1/(2*epsilon)
    epsq2 = 1/(2*epsilon^2)
    vshift = 2*cellnumber
    B = zeros(size(B))
    if fluxtype == "full"
        for i in 1:2*cellnumber
            # auf rho wirkende Terme
            if i <= cellnumber
                # orange/lila
                if i != 1
                    B[2*i-1, 2*(i-1)] = -eps2
                    B[2*i-1, 2*(i-1)+vshift] = -1/2
                end
                B[1, 2*cellnumber] = -eps2
                B[1, 2*cellnumber+vshift] = -1/2
                #grün/gelb
                B[2*i-1, 2*i-1] = eps2
                B[2*i-1, 2*i-1+vshift] = -1/2
                #rot/pink
                B[2*i, 2*i] = eps2
                B[2*i, 2*i+vshift] = 1/2
                #blau/braun
                if i != cellnumber
                    B[2*i, 2*i+1] = -eps2
                    B[2*i, 2*i+1+vshift] = 1/2
                end
                B[2*cellnumber, 1] = -eps2
                B[2*cellnumber, 1+vshift] = 1/2
            # auf j wirkende Terme
            else
                # orange/lila
                if i != 1 + cellnumber
                    B[2*i-1, 2*(i-1)-vshift] = -epsq2
                    B[2*i-1, 2*(i-1)] = -eps2
                end
                B[1+vshift, 2*cellnumber] = -epsq2
                B[1+vshift, 2*cellnumber+vshift] = -eps2
                #grün/gelb
                B[2*i-1, 2*i-1-vshift] = -epsq2
                B[2*i-1, 2*i-1] = eps2
                #rot/pink
                B[2*i, 2*i-vshift] = epsq2
                B[2*i, 2*i] = eps2
                #blau/braun
                if i != 2*cellnumber
                    B[2*i, 2*i+1-vshift] = epsq2
                    B[2*i, 2*i+1] = -eps2
                end
                B[2*cellnumber+vshift, 1] = epsq2
                B[2*cellnumber+vshift, 1+vshift] = -eps2
            end
        end
        elseif fluxtype == "full2"
            for i in 1:2*cellnumber
                # auf rho wirkende Terme
                if i <= cellnumber
                    # orange/lila
                    if i != 1
                        B[2*i-1, 2*(i-1)] = -eps2
                        B[2*i-1, 2*(i-1)+vshift] = -epsq2
                    end
                    B[1, 2*cellnumber] = -eps2
                    B[1, 2*cellnumber+vshift] = -epsq2
                    #grün/gelb
                    B[2*i-1, 2*i-1] = eps2
                    B[2*i-1, 2*i-1+vshift] = -epsq2
                    #rot/pink
                    B[2*i, 2*i] = eps2
                    B[2*i, 2*i+vshift] = epsq2
                    #blau/braun
                    if i != cellnumber
                        B[2*i, 2*i+1] = -eps2
                        B[2*i, 2*i+1+vshift] = epsq2
                    end
                    B[2*cellnumber, 1] = -eps2
                    B[2*cellnumber, 1+vshift] = epsq2
                # auf j wirkende Terme
                else
                    # orange/lila
                    if i != 1 + cellnumber
                        B[2*i-1, 2*(i-1)-vshift] = -1/2
                        B[2*i-1, 2*(i-1)] = -eps2
                    end
                    B[1+vshift, 2*cellnumber] = -1/2
                    B[1+vshift, 2*cellnumber+vshift] = -eps2
                    #grün/gelb
                    B[2*i-1, 2*i-1-vshift] = -1/2
                    B[2*i-1, 2*i-1] = eps2
                    #rot/pink
                    B[2*i, 2*i-vshift] = 1/2
                    B[2*i, 2*i] = eps2
                    #blau/braun
                    if i != 2*cellnumber
                        B[2*i, 2*i+1-vshift] = 1/2
                        B[2*i, 2*i+1] = -eps2
                    end
                    B[2*cellnumber+vshift, 1] = 1/2
                    B[2*cellnumber+vshift, 1+vshift] = -eps2
                end
            end
    elseif fluxtype == "altlr"
        for i in 1:2*cellnumber
            # auf rho wirkende Terme
            if i <= cellnumber
                # lila
                if i != 1
                    B[2*i-1, 2*(i-1)+vshift] = -1/2 *2
                end
                B[1, 2*cellnumber+vshift] = -1/2 *2
                # pink
                B[2*i, 2*i+vshift] = 1/2 *2
            # auf j wirkende Terme
            else
                # grün
                B[2*i-1, 2*i-1-vshift] = -epsq2 *2
                # blau
                if i != 2*cellnumber
                    B[2*i, 2*i+1-vshift] = epsq2 *2
                end
                B[2*cellnumber+vshift, 1] = epsq2 *2
            end
        end
    elseif fluxtype == "altrl"
        for i in 1:2*cellnumber
            # auf rho wirkende Terme
            if i <= cellnumber
                #gelb
                B[2*i-1, 2*i-1+vshift] = -1/2 *2
                #braun
                if i != cellnumber
                    B[2*i, 2*i+1+vshift] = 1/2 *2
                end
                B[2*cellnumber, 1+vshift] = 1/2 *2
            # auf j wirkende Terme
            else
                # orange/lila
                if i != 1 + cellnumber
                    B[2*i-1, 2*(i-1)-vshift] = -epsq2 *2
                end
                B[1+vshift, 2*cellnumber] = -epsq2 *2
                #rot
                B[2*i, 2*i-vshift] = epsq2 *2
            end
        end
    else
        println("No numerical flux chosen! Therefore, B=0")
    end
    return B
end
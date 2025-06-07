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
        elseif init_cond == "eq_data"
            return sin(x)
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
        if eq_type == "wave"
            u_wave_exact = zeros((deg+1)*2*N, length(t_d));
            u_wave_exact[1:(deg+1)*N, :] = [1/2*(sin(2*pi*(x-1/epsilon*t))+epsilon*cos(2*pi*(x-1/epsilon*t))+sin(2*pi*(x+1/epsilon*t))-epsilon*cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in t_d];
            u_wave_exact[(deg+1)*N+1:(deg+1)*2*N, :] = [1/2*(1/epsilon*sin(2*pi*(x-1/epsilon*t))+cos(2*pi*(x-1/epsilon*t))-1/epsilon*sin(2*pi*(x+1/epsilon*t))+cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in t_d];
            return u_wave_exact
        elseif eq_type == "telegraph"
            if setup["init_cond"] != "telsin"
                println("Care: exact solution does not match inital conditions!")
            end
            r = -2/(1+sqrt(1-4*epsilon^2))
            u_tel_exact = zeros((deg+1)*2*N, length(t_d));
            u_tel_exact[1:(deg+1)*N, :] = [1/r*exp(r*t)*sin(x) for x in x_d, t in t_d];
            u_tel_exact[(deg+1)*N+1:(deg+1)*2*N, :] = [exp(r*t)*cos(x) for x in x_d, t in t_d];
            return u_tel_exact
        elseif eq_type == "heat"
            r = -1
            u_tel_exact = zeros((deg+1)*N, length(t_d));
            u_tel_exact[1:(deg+1)*N, :] = [1/r*exp(r*t)*sin(x) for x in x_d, t in t_d];
        else
            println("wrong eq type")
        end
    end
end

function f_RHS(u, t, RHS_mat, setup)
    #u0_eval = setup["u0_eval"]
    #x_d = setup["x_d"]
    #a = setup["a"]
    #basis = setup["basis"]
    #deg = setup["deg"]
    #v = setup["v"]
    # Nx = length(x_d)*2, because wave/telegraph is a 2-dim system of equations
    #Nx = length(x_d)*2
    return RHS_mat * u
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

function DGsemidiscretization_DoD_telegraph(setup, deg, nodes, num_flux;
                                             T = Float64, do_stabilize = false, c = 1.0, fix_eta = true, eq_type = "telegraph", fluxtype = "full", include_b_h = true, ext_test_func = true)
    cellnumber = setup["cellnumber"]
    v = setup["v"]
    a = setup["a"]
    u0_eval = setup["u0_eval"]
    t_d = setup["t_d"]
    h = setup["h"]
    init_cond = setup["init_cond"]
    setup["nodes"] = nodes
    setup["deg"] = deg
    setup["num_flux"] = num_flux
    setup["fix_eta"] = fix_eta
    setup["c"] = c

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
    setup["xrange"] = xrange
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

    ex_DM = zeros(T, xrange, xrange)
    ex_sM = zeros(T, xrange, xrange)
    ex_B = zeros(T, cellnumber*4, cellnumber*4)

    im_DM = zeros(T, xrange, xrange)
    im_sM = zeros(T, xrange, xrange)
    im_B = zeros(T, cellnumber*4, cellnumber*4)

    A_DM = zeros(T, xrange, xrange)
    A_B = zeros(T, cellnumber*4, cellnumber*4)

    ex_DM_J1 = zeros(T, xrange, xrange)
    im_DM_J1 = zeros(T, xrange, xrange)
    im_DM_J1 = zeros(T, xrange, xrange)    
    ex_B_J0 = zeros(T, cellnumber*4, cellnumber*4)
    ex_B_J0E1 = zeros(T, cellnumber*4, cellnumber*4)
    ex_B_J0E2 = zeros(T, cellnumber*4, cellnumber*4)
    im_B_J0 = zeros(T, cellnumber*4, cellnumber*4)
    im_B_J0E1 = zeros(T, cellnumber*4, cellnumber*4)
    im_B_J0E2 = zeros(T, cellnumber*4, cellnumber*4)
    A_B_J0 = zeros(T, cellnumber*4, cellnumber*4)
    A_B_J0E1 = zeros(T, cellnumber*4, cellnumber*4)
    A_B_J0E2 = zeros(T, cellnumber*4, cellnumber*4)
    V = zeros(T, cellnumber*4, xrange)
    ex_V_J0E1 = zeros(T, cellnumber*4, xrange)
    ex_V_J0E2 = zeros(T, cellnumber*4, xrange)
    im_V_J0E1 = zeros(T, cellnumber*4, xrange)
    im_V_J0E2 = zeros(T, cellnumber*4, xrange)
    A_V_J0E1 = zeros(T, cellnumber*4, xrange)
    A_V_J0E2 = zeros(T, cellnumber*4, xrange)

    ex_A = zeros(T, xrange, xrange)

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
            ex_DM[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(cellnumber+(i-1)) + 1:(deg+1)*(cellnumber+i)] = basis.D' * diagm(basis.weights)
        else
            im_DM[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*((i-1)-cellnumber) + 1:(deg+1)*(i-cellnumber)] = basis.D' * diagm(basis.weights)/(epsilon^2)
            im_sM[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] = diagm(basis.weights)/(epsilon^2)
        end

        # Constructing the Interpolation to the vertices for the flux terms with Gauss-Legendre nodes
        #V[2*i-1, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i], basis1)
        #V[2*i, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+1], basis1)
        V[2*i-1, (deg+1)*(i-1) + 1:(deg+1)*i] = left_interpolation
        V[2*i, (deg+1)*(i-1) + 1:(deg+1)*i] = right_interpolation
        
        # Extrapolation to the Inflow of Cell + 1(V_J0E1) and the Outflow of Cut-Cell(V_J0E2) from Cell - 1 (from the perspective of Cell -1).
        #=
        if do_stabilize
            if i < cellnumber
                if 2*(v[i+2]-v[i+1])/h<1
                    # Dies ist nur für den altlr Fluss implementiert. Für Fluss in andere Richtung, müsste hier eine Fallunterscheidung einbaut werden
                    # (und selbstverständlich die B_J1/J2 Terme angepasst werden). Andere mögliche Einträge für V_J0E1/E2 sind aktuell auskommentiert
                    #V_J0E1[2*i+2, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+2], basis1) 
                    #V_J0E1[2*i+2 + 2*cellnumber, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+2], basis1) 
                    #V_J0E1[2*i+2, (deg+1)*(i-1) + 1 + cellnumber:(deg+1)*i + cellnumber] = interpolation_matrix(v[i+2], basis1) 
                    ex_V_J0E1[2*i+2 + 2*cellnumber, (deg+1)*(i-1) + 1 + cellnumber:(deg+1)*i + cellnumber] = interpolation_matrix(v[i+2], basis1) 

                    #V_J0E2[2*i+3, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+2], basis1)
                    #V_J0E2[2*i+3 + 2*cellnumber, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+2], basis1)
                    #V_J0E2[2*i+3, (deg+1)*(i-1) + 1 + cellnumber:(deg+1)*i + cellnumber] = interpolation_matrix(v[i+2], basis1)
                    ex_V_J0E2[2*i+3 + 2*cellnumber, (deg+1)*(i-1) + 1 + cellnumber:(deg+1)*i + cellnumber] = interpolation_matrix(v[i+2], basis1)

                    #im_V_J0E1[2*i+2 + 2*cellnumber, (deg+1)*(i-1) + 1 + cellnumber:(deg+1)*i + cellnumber] = interpolation_matrix(v[i+2], basis1) 
                    #im_V_J0E1[2*i+2 + 2*cellnumber, (deg+1)*(i-1) + 1 + cellnumber:(deg+1)*i + cellnumber] = interpolation_matrix(v[i+2], basis1) 
                end
            end
        end
        =#
        #
        if do_stabilize
            if i < cellnumber
                if 2*(v[i+1]-v[i])/h<1
                    if fluxtype == "altlr"
                         # cell - 1 is inflow cell, therefore construct basis1 on cell - 1
                        for j in 1:deg+1
                            basis1.nodes[j] = x_d[(deg+1)*(i-2) + j] 
                        end
                        # Dies ist nur für den altlr Fluss implementiert. Für Fluss in andere Richtung, müsste hier eine Fallunterscheidung einbaut werden
                        # (und selbstverständlich die B_J1/J2 Terme angepasst werden). Andere mögliche Einträge für V_J0E1/E2 sind aktuell auskommentiert
                        #extrapolation to right flow
                        ex_V_J0E1[2*(i-1)+2 + 2*cellnumber, (deg+1)*(cellnumber + (i-2)) + 1:(deg+1)*(cellnumber + (i-1))] = interpolation_matrix(v[i+1], basis1) 
                        ex_V_J0E2[2*(i-1)+3 + 2*cellnumber, (deg+1)*(cellnumber + (i-2)) + 1:(deg+1)*(cellnumber + (i-1))] = interpolation_matrix(v[i+1], basis1)
                        # cell + 1 is inflow cell, therefore construct basis1 on cell + 1
                        for j in 1:deg+1
                            basis1.nodes[j] = x_d[(deg+1)*(i) + j] 
                        end
                        im_V_J0E1[2*i-1, (deg+1)*i + 1:(deg+1)*(i+1)] = interpolation_matrix(v[i], basis1) 
                        im_V_J0E2[2*i-2, (deg+1)*i + 1:(deg+1)*(i+1)] = interpolation_matrix(v[i], basis1)
                    elseif fluxtype == "full"
                        for j in 1:deg+1
                            basis1.nodes[j] = x_d[(deg+1)*(i-2) + j] 
                        end
                        #bzgl explizit, pink
                        ex_V_J0E1[2*(i-1)+2 + 2*cellnumber, (deg+1)*(cellnumber + (i-2)) + 1:(deg+1)*(cellnumber + (i-1))] = interpolation_matrix(v[i+1], basis1) 
                        ex_V_J0E2[2*(i-1)+3 + 2*cellnumber, (deg+1)*(cellnumber + (i-2)) + 1:(deg+1)*(cellnumber + (i-1))] = interpolation_matrix(v[i+1], basis1)
                        #bzgl explizit, rot
                        ex_V_J0E1[2*(i-1)+2, (deg+1)*(i-2) + 1:(deg+1)*(i-1)] = interpolation_matrix(v[i+1], basis1) 
                        ex_V_J0E2[2*(i-1)+3, (deg+1)*(i-2) + 1:(deg+1)*(i-1)] = interpolation_matrix(v[i+1], basis1)
                        #bzgl implizit, pink
                        im_V_J0E1[2*(i-1)+2 + 2*cellnumber, (deg+1)*(cellnumber + (i-2)) + 1:(deg+1)*(cellnumber + (i-1))] = interpolation_matrix(v[i+1], basis1) 
                        im_V_J0E2[2*(i-1)+3 + 2*cellnumber, (deg+1)*(cellnumber + (i-2)) + 1:(deg+1)*(cellnumber + (i-1))] = interpolation_matrix(v[i+1], basis1)
                        #bzgl implizit, rot
                        im_V_J0E1[2*(i-1)+2, (deg+1)*(i-2) + 1:(deg+1)*(i-1)] = interpolation_matrix(v[i+1], basis1) 
                        im_V_J0E2[2*(i-1)+3, (deg+1)*(i-2) + 1:(deg+1)*(i-1)] = interpolation_matrix(v[i+1], basis1)
                        for j in 1:deg+1
                            basis1.nodes[j] = x_d[(deg+1)*(i) + j] 
                        end
                        #bzgl. explizit, gelb
                        ex_V_J0E1[2*i-1 + 2*cellnumber, (deg+1)*(cellnumber + (i)) + 1:(deg+1)*(cellnumber + (i+1))] = interpolation_matrix(v[i], basis1)
                        ex_V_J0E2[2*i-2 + 2*cellnumber, (deg+1)*(cellnumber + (i)) + 1:(deg+1)*(cellnumber + (i+1))] = interpolation_matrix(v[i], basis1)
                        #bzgl. explizit, grün
                        ex_V_J0E1[2*i-1, (deg+1)*i + 1:(deg+1)*(i+1)] = interpolation_matrix(v[i], basis1) 
                        ex_V_J0E2[2*i-2, (deg+1)*i + 1:(deg+1)*(i+1)] = interpolation_matrix(v[i], basis1)
                        #bzgl. implizit, gelb
                        im_V_J0E1[2*i-1 + 2*cellnumber, (deg+1)*(cellnumber + (i)) + 1:(deg+1)*(cellnumber + (i+1))] = interpolation_matrix(v[i], basis1)
                        im_V_J0E2[2*i-2 + 2*cellnumber, (deg+1)*(cellnumber + (i)) + 1:(deg+1)*(cellnumber + (i+1))] = interpolation_matrix(v[i], basis1)
                        #bzgl. implizit, grün
                        im_V_J0E1[2*i-1, (deg+1)*i + 1:(deg+1)*(i+1)] = interpolation_matrix(v[i], basis1) 
                        im_V_J0E2[2*i-2, (deg+1)*i + 1:(deg+1)*(i+1)] = interpolation_matrix(v[i], basis1)
                    end
                    # ACHTUNG!! Neben den im_VJ0E1/im_VJ0E2 termen ist es auch wichtig nochmal auf den Term zu schauen, der uns gerade die Domain-of-Dependence in der 2. Komponente erweitert,
                    # ohne das wir es möchten! Eigentlich sollten an dieser stelle ja noch keine extra Terme in der semidiskretisierung auftreten, sondern nur terme angepasst werden
                end
            end
        end
        #
            
        TranspV[(deg+1)*(i-1)+1:(deg+1)*i, 2*i-1:2*i] = V[2*i-1:2*i, (deg+1)*(i-1) + 1:(deg+1)*i]'
    end

    # Constructing numerical fluxes(with periodic BC)
    ex_B, im_B = flux_terms_telegraph!(ex_B, im_B, setup, deg, fluxtype = fluxtype)
    # DoD-Stabilization
    if do_stabilize
        ex_B_J0, ex_B_J0E1, ex_B_J0E2, im_B_J0, im_B_J0E1, im_B_J0E2 = flux_terms_telegraph_DoD!(ex_B, im_B, setup, deg, fluxtype = fluxtype)
    end

    # DoD-Stabilization
    if do_stabilize
        for i in 1:cellnumber
            if 2*(v[i+1]-v[i])/h<1
                if fix_eta == false
                    eta = 1-minimum([1/(2*deg+1)*(v[i+1]-v[i])/((t_d[2]-t_d[1])*a), 1.0])
                else
                    eta = 1-minimum([(v[i+1]-v[i])/h*1/c, 1.0])
                end
                # Volume Term
                if deg>0
                    ex_DM_J1[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(cellnumber + (i-1)) + 1:(deg+1)*(cellnumber + i)] +=
                        - basis.D' * diagm(basis.weights) * eta
                    
                    basis2=nodes(deg)
                    for j in 1:deg+1
                        # Inserting these nodes into the structure basis2, to adapt the field basis1.interp_matrix to the used nodes.
                        basis2.nodes[j] = x_d[(deg+1)*(i-2) + j] 
                    end
                    Vol_interpol = interpolation_matrix(x_d[(deg+1)*(i-1) + 1:(deg+1)*(i)] , basis2)
                    ex_DM_J1[(deg+1)*(i-1) + 1:(deg+1)*(i), (deg+1)*(cellnumber + (i-2)) + 1:(deg+1)*(cellnumber + (i-1))] += 
                        + basis.D' * diagm(basis.weights) * Vol_interpol * eta
                    # Extra-Stabilityterm with the extended test-function (L_E_in(u_h)-u_h)*L_E_in(v_h)
                    if ext_test_func == true
                        ex_DM_J1[(deg+1)*(i-2) + 1:(deg+1)*(i-1), (deg+1)*(cellnumber + (i-1)) + 1:(deg+1)*(cellnumber + i)] +=
                            Vol_interpol' *  basis.D' * diagm(basis.weights) * eta
                        ex_DM_J1[(deg+1)*(i-2) + 1:(deg+1)*(i-1), (deg+1)*(cellnumber + (i-2)) + 1:(deg+1)*(cellnumber + (i-1))] += 
                            - Vol_interpol' * basis.D' * diagm(basis.weights) * Vol_interpol * eta
                    end


                    im_DM_J1[(deg+1)*(cellnumber + (i-1)) + 1:(deg+1)*(cellnumber + i), (deg+1)*(i-1) + 1:(deg+1)*i] +=
                        - basis.D' * diagm(basis.weights) * eta /(epsilon^2)
                    
                    basis2=nodes(deg)
                    for j in 1:deg+1
                        # Inserting these nodes into the structure basis2, to adapt the field basis1.interp_matrix to the used nodes.
                        basis2.nodes[j] = x_d[(deg+1)*(i) + j] 
                    end
                    Vol_interpol = interpolation_matrix(x_d[(deg+1)*(i-1) + 1:(deg+1)*(i)] , basis2)
                    im_DM_J1[(deg+1)*(cellnumber + (i-1)) + 1:(deg+1)*(cellnumber + i), (deg+1)*(i) + 1:(deg+1)*(i+1)] += 
                        + basis.D' * diagm(basis.weights) * Vol_interpol * eta /(epsilon^2)
                    # Extra-Stabilityterm with the extended test-function (L_E_in(u_h)-u_h)*L_E_in(v_h)
                    if ext_test_func == true
                        im_DM_J1[(deg+1)*(cellnumber + (i)) + 1:(deg+1)*(cellnumber + (i+1)), (deg+1)*(i-1) + 1:(deg+1)*i] +=
                            Vol_interpol' *  basis.D' * diagm(basis.weights) * eta /(epsilon^2)
                        im_DM_J1[(deg+1)*(cellnumber + i) + 1:(deg+1)*(cellnumber + (i+1)), (deg+1)*(i) + 1:(deg+1)*(i+1)] += 
                            - Vol_interpol' * basis.D' * diagm(basis.weights) * Vol_interpol * eta /(epsilon^2)
                    end
                end
            end
        end
    end
    
    # Constructing RHS Matrix
    invM = inv(M)
    Mh = inv(Sc[1:Int(xrange/2), 1:Int(xrange/2)])*M[1:Int(xrange/2), 1:Int(xrange/2)]
    if eq_type == "telegraph"
        splitRHS["Dminus"] = -(Sc*invM*(-TranspV* (ex_B*V + ex_B_J0*V + ex_B_J0E1*ex_V_J0E1 + ex_B_J0E2*ex_V_J0E2) + ex_DM + ex_DM_J1))[1:Int(xrange/2), Int(xrange/2)+1:xrange]
        splitRHS["Dplus"] = -epsilon^2*(Sc*invM*(-TranspV* (im_B*V + im_B_J0*V + im_B_J0E1*im_V_J0E1 + im_B_J0E2*im_V_J0E2) + im_DM + im_DM_J1))[Int(xrange/2)+1:xrange, 1:Int(xrange/2)]
        splitRHS["Dminus2"] = -(Sc*invM*(-TranspV* (ex_B*V + ex_B_J0*V + ex_B_J0E1*ex_V_J0E1 + ex_B_J0E2*ex_V_J0E2) + ex_DM ))[1:Int(xrange/2), Int(xrange/2)+1:xrange]
        splitRHS["Dplus2"] = -epsilon^2*(Sc*invM*(-TranspV* (im_B*V + im_B_J0*V + im_B_J0E1*im_V_J0E1 + im_B_J0E2*im_V_J0E2) + im_DM))[Int(xrange/2)+1:xrange, 1:Int(xrange/2)]
        A = 1/(2*epsilon)*(splitRHS["Dplus"]-splitRHS["Dminus"])
        splitRHS["A"] = A
        if include_b_h == true
            ex_A[Int(xrange/2)+1:xrange, Int(xrange/2)+1:xrange] = A
        end
        splitRHS["IntByParts"] = Mh*splitRHS["Dplus"] + splitRHS["Dminus"]'*Mh
        splitRHS["symm_negsemi"] = Mh*(splitRHS["Dplus"]-splitRHS["Dminus"])

        ex_RHS_mat = Sc*invM * (-TranspV* (ex_B*V + ex_B_J0*V + ex_B_J0E1*ex_V_J0E1 + ex_B_J0E2*ex_V_J0E2) + ex_DM + ex_DM_J1) -invM*ex_sM + ex_A
        im_RHS_mat = Sc*invM * (-TranspV* (im_B*V + im_B_J0*V + im_B_J0E1*im_V_J0E1 + im_B_J0E2*im_V_J0E2) + im_DM + im_DM_J1) -invM*im_sM
        RHS_mat = Sc*invM * (-TranspV* ((ex_B + im_B)*V + (ex_B_J0+im_B_J0)*V + ex_B_J0E1*ex_V_J0E1 + im_B_J0E1*im_V_J0E1 + ex_B_J0E2*ex_V_J0E2 + im_B_J0E2*im_V_J0E2) + (ex_DM+im_DM) + (ex_DM_J1 + im_DM_J1)) -invM*(ex_sM + im_sM) + ex_A
    elseif eq_type == "wave"
        #Wave equation (choosing epsilon =1/c), not applicable for DoD here
        splitRHS["A"] = zeros(xrange, xrange)
        ex_RHS_mat = Sc*invM * (-TranspV* (ex_B*V + ex_B_J0*V + ex_B_J0E1*ex_V_J0E1 + ex_B_J0E2*ex_V_J0E2) + ex_DM + ex_DM_J1)
        im_RHS_mat = Sc*invM * (-TranspV* (im_B*V + im_B_J0*V + im_B_J0E1*im_V_J0E1 + im_B_J0E2*im_V_J0E2) + im_DM + im_DM_J1)
        RHS_mat = Sc*invM * (-TranspV* ((ex_B + im_B)*V + (ex_B_J0+im_B_J0)*V + ex_B_J0E1*ex_V_J0E1 + im_B_J0E1*im_V_J0E1 + ex_B_J0E2*ex_V_J0E2 + im_B_J0E2*im_V_J0E2) + (ex_DM+im_DM) + (ex_DM_J1 + im_DM_J1))
    elseif eq_type == "heat"
        if fluxtype == "altlr"
            RHS_mat = DGsemidiscretization_DoD_heat(setup, deg, nodes, num_flux; T = T, do_stabilize = do_stabilize, c = c, fix_eta = fix_eta, ext_test_func = ext_test_func)[1]
            splitRHS["Dminus"] = -(Sc*invM*(-TranspV* (ex_B*V + ex_B_J0*V + ex_B_J0E1*ex_V_J0E1 + ex_B_J0E2*ex_V_J0E2) + ex_DM + ex_DM_J1))[1:Int(xrange/2), Int(xrange/2)+1:xrange]
            splitRHS["Dplus"] = -epsilon^2*(Sc*invM*(-TranspV* (im_B*V + im_B_J0*V + im_B_J0E1*im_V_J0E1 + im_B_J0E2*im_V_J0E2) + im_DM + im_DM_J1))[Int(xrange/2)+1:xrange, 1:Int(xrange/2)]
            ## evtl einfacher???:
            #RHS_mat = splitRHS["Dminus"]*splitRHS["Dplus"] (muss aus eq_type == "telegraph" gezogen werden -> einfach immer aufrufen und am ende ex_RHS_mat, im_RHS_mat, RHS_mat überschreiben?)
            ex_RHS_mat = RHS_mat
            im_RHS_mat = zeros(size(RHS_mat))
        else
            throw(ArgumentError("That flux is not implemented for the heat equation"))
        end
    end

    setup["ex_RHS_mat"] = ex_RHS_mat
    setup["im_RHS_mat"] = im_RHS_mat
    setup["RHS_mat"] = RHS_mat

    # Store parts in dictionary
    splitRHS["Sc"] = Sc
    splitRHS["M"] = M
    splitRHS["invM"] = invM
    splitRHS["TranspV"] = TranspV
    splitRHS["V"] = V
    
    splitRHS["ex_sM"] = ex_sM
    splitRHS["im_sM"] = im_sM
    splitRHS["ex_B"] = ex_B
    splitRHS["im_B"] = im_B
    splitRHS["ex_DM"] = ex_DM
    splitRHS["im_DM"] = im_DM
    splitRHS["TranspV*ex_B*V"] = TranspV*ex_B*V
    splitRHS["TranspV*im_B*V"] = TranspV*im_B*V

    splitRHS["ex_V_J0E1"] = ex_V_J0E1
    splitRHS["ex_V_J0E2"] = ex_V_J0E2
    splitRHS["im_V_J0E1"] = im_V_J0E1
    splitRHS["im_V_J0E2"] = im_V_J0E2
    splitRHS["ex_B_J0"] = ex_B_J0
    splitRHS["ex_B_J0E1"] = ex_B_J0E1
    splitRHS["ex_B_J0E2"] = ex_B_J0E2
    splitRHS["im_B_J0"] = im_B_J0
    splitRHS["im_B_J0E1"] = im_B_J0E1
    splitRHS["im_B_J0E2"] = im_B_J0E2
    splitRHS["TranspV*ex_B_J0*V"] = TranspV*ex_B_J0*V
    splitRHS["TranspV*im_B_J0*V"] = TranspV*im_B_J0*V
    splitRHS["TranspV*ex_B_J0E1*ex_V_J0E1"] = TranspV*ex_B_J0E1*ex_V_J0E1
    splitRHS["TranspV*im_B_J0E1*im_V_J0E1"] = TranspV*im_B_J0E1*im_V_J0E1
    splitRHS["TranspV*ex_B_J0E2*ex_V_J0E2"] = TranspV*ex_B_J0E2*ex_V_J0E2
    splitRHS["TranspV*im_B_J0E2*im_V_J0E2"] = TranspV*im_B_J0E2*im_V_J0E2
    splitRHS["ex_DM_J1"] = ex_DM_J1
    splitRHS["im_DM_J1"] = im_DM_J1

    splitRHS["FLUXTERMS"] = Sc*invM * (-TranspV* ((ex_B + im_B)*V + (ex_B_J0+im_B_J0)*V + ex_B_J0E1*ex_V_J0E1 + im_B_J0E1*im_V_J0E1 + ex_B_J0E2*ex_V_J0E2 + im_B_J0E2*im_V_J0E2))
    splitRHS["VOLTERMS"] = Sc*invM * ((ex_DM+im_DM) + (ex_DM_J1 + im_DM_J1)) -invM*(ex_sM + im_sM)
    splitRHS["Mh"] = Mh
    # Constructing the initial conition (special case for j coordinate)
    u0_rho = [u0_eval(x) for x in x_d]
    if eq_type == "heat"
        if init_cond == "telsin"
            u0 = [-sin(x) for x in x_d]
        elseif init_cond == "eq_data"
            u0 = [sin(x) for x in x_d]
        end
        M_global = M_global[1:Int(xrange/2), 1:Int(xrange/2)]
    else
        if eq_type == "telegraph"
            if init_cond == "telsin"
                u0_j = cos.(x_d)
            elseif init_cond == "eq_data"
                u0_j = cos.(x_d)
            end
        elseif eq_type == "wave"
            u0_j = cos.(2*pi*x_d)
        end
        u0 = vcat(u0_rho, u0_j);
    end


    setup["u0"] = u0
    setup["M"] = M_global
    setup["x_d"] = x_d
    setup["basis"] = basis
    #display((-TranspV* (B*V + B2*V2))[1:9, 1:9])
    #display(DM_J1[1:10, 1:10])
    
    return RHS_mat, setup, splitRHS
end


function flux_terms_telegraph!(ex_B, im_B, setup, deg; fluxtype = "altlr")
    cellnumber = setup["cellnumber"]
    epsilon = setup["epsilon"]
    eps2 = 1/(2*epsilon)
    epsq2 = 1/(2*epsilon^2)
    vshift = 2*cellnumber
    ex_B = zeros(size(ex_B))
    im_B = zeros(size(im_B))
    if fluxtype == "full"
        for i in 1:2*cellnumber
            # auf rho wirkende Terme
            if i <= cellnumber
                # orange/lila
                if i != 1
                    ex_B[2*i-1, 2*(i-1)] = -eps2
                    ex_B[2*i-1, 2*(i-1)+vshift] = -1/2
                end
                ex_B[1, 2*cellnumber] = -eps2
                ex_B[1, 2*cellnumber+vshift] = -1/2
                #grün/gelb
                ex_B[2*i-1, 2*i-1] = eps2
                ex_B[2*i-1, 2*i-1+vshift] = -1/2
                #rot/pink
                ex_B[2*i, 2*i] = eps2
                ex_B[2*i, 2*i+vshift] = 1/2
                #blau/braun
                if i != cellnumber
                    ex_B[2*i, 2*i+1] = -eps2
                    ex_B[2*i, 2*i+1+vshift] = 1/2
                end
                ex_B[2*cellnumber, 1] = -eps2
                ex_B[2*cellnumber, 1+vshift] = 1/2
            # auf j wirkende Terme
            else
                # orange/lila
                if i != 1 + cellnumber
                    im_B[2*i-1, 2*(i-1)-vshift] = -epsq2
                    im_B[2*i-1, 2*(i-1)] = -eps2
                end
                im_B[1+vshift, 2*cellnumber] = -epsq2
                im_B[1+vshift, 2*cellnumber+vshift] = -eps2
                #grün/gelb
                im_B[2*i-1, 2*i-1-vshift] = -epsq2
                im_B[2*i-1, 2*i-1] = eps2
                #rot/pink
                im_B[2*i, 2*i-vshift] = epsq2
                im_B[2*i, 2*i] = eps2
                #blau/braun
                if i != 2*cellnumber
                    im_B[2*i, 2*i+1-vshift] = epsq2
                    im_B[2*i, 2*i+1] = -eps2
                end
                im_B[2*cellnumber+vshift, 1] = epsq2
                im_B[2*cellnumber+vshift, 1+vshift] = -eps2
            end
        end
        elseif fluxtype == "full2"
            for i in 1:2*cellnumber
                # auf rho wirkende Terme
                if i <= cellnumber
                    # orange/lila
                    if i != 1
                        ex_B[2*i-1, 2*(i-1)] = -eps2
                        ex_B[2*i-1, 2*(i-1)+vshift] = -epsq2
                    end
                    ex_B[1, 2*cellnumber] = -eps2
                    ex_B[1, 2*cellnumber+vshift] = -epsq2
                    #grün/gelb
                    ex_B[2*i-1, 2*i-1] = eps2
                    ex_B[2*i-1, 2*i-1+vshift] = -epsq2
                    #rot/pink
                    ex_B[2*i, 2*i] = eps2
                    ex_B[2*i, 2*i+vshift] = epsq2
                    #blau/braun
                    if i != cellnumber
                        ex_B[2*i, 2*i+1] = -eps2
                        ex_B[2*i, 2*i+1+vshift] = epsq2
                    end
                    ex_B[2*cellnumber, 1] = -eps2
                    ex_B[2*cellnumber, 1+vshift] = epsq2
                # auf j wirkende Terme
                else
                    # orange/lila
                    if i != 1 + cellnumber
                        im_B[2*i-1, 2*(i-1)-vshift] = -1/2
                        im_B[2*i-1, 2*(i-1)] = -eps2
                    end
                    im_B[1+vshift, 2*cellnumber] = -1/2
                    im_B[1+vshift, 2*cellnumber+vshift] = -eps2
                    #grün/gelb
                    im_B[2*i-1, 2*i-1-vshift] = -1/2
                    im_B[2*i-1, 2*i-1] = eps2
                    #rot/pink
                    im_B[2*i, 2*i-vshift] = 1/2
                    im_B[2*i, 2*i] = eps2
                    #blau/braun
                    if i != 2*cellnumber
                        im_B[2*i, 2*i+1-vshift] = 1/2
                        im_B[2*i, 2*i+1] = -eps2
                    end
                    im_B[2*cellnumber+vshift, 1] = 1/2
                    im_B[2*cellnumber+vshift, 1+vshift] = -eps2
                end
            end
    elseif fluxtype == "altlr"
        for i in 1:2*cellnumber
            # auf rho wirkende Terme
            if i <= cellnumber
                # lila
                if i != 1
                    ex_B[2*i-1, 2*(i-1)+vshift] = -1/2 *2
                end
                ex_B[1, 2*cellnumber+vshift] = -1/2 *2
                # pink
                ex_B[2*i, 2*i+vshift] = 1/2 *2
            # auf j wirkende Terme
            else
                # grün
                im_B[2*i-1, 2*i-1-vshift] = -epsq2 *2
                # blau
                if i != 2*cellnumber
                    im_B[2*i, 2*i+1-vshift] = epsq2 *2
                end
                im_B[2*cellnumber+vshift, 1] = epsq2 *2
            end
        end
    elseif fluxtype == "altrl"
        for i in 1:2*cellnumber
            # auf rho wirkende Terme
            if i <= cellnumber
                #gelb
                ex_B[2*i-1, 2*i-1+vshift] = -1/2 *2
                #braun
                if i != cellnumber
                    ex_B[2*i, 2*i+1+vshift] = 1/2 *2
                end
                ex_B[2*cellnumber, 1+vshift] = 1/2 *2
            # auf j wirkende Terme
            else
                # orange/lila
                if i != 1 + cellnumber
                    im_B[2*i-1, 2*(i-1)-vshift] = -epsq2 *2
                end
                im_B[1+vshift, 2*cellnumber] = -epsq2 *2
                #rot
                im_B[2*i, 2*i-vshift] = epsq2 *2
            end
        end
    else
        println("No numerical flux chosen! Therefore, B=0")
    end
    return ex_B, im_B
end


function flux_terms_telegraph_DoD!(ex_B, im_B, setup, deg; fluxtype = "altlr")
    cellnumber = setup["cellnumber"]
    epsilon = setup["epsilon"]
    eps2 = 1/(2*epsilon)
    epsq2 = 1/(2*epsilon^2)
    v = setup["v"]
    a = setup["a"]
    t_d = setup["t_d"]
    h = setup["h"]
    fix_eta = setup["fix_eta"]
    c = setup["c"]
    vshift = 2*cellnumber
    ex_B_J0 = zeros(size(ex_B))
    ex_B_J0E1 = zeros(size(ex_B))
    ex_B_J0E2 = zeros(size(ex_B))
    im_B_J0 = zeros(size(im_B))
    im_B_J0E1 = zeros(size(im_B))
    im_B_J0E2 = zeros(size(im_B))
    for i in 1:cellnumber
        if 2*(v[i+1]-v[i])/h<1
            if fix_eta == false
                eta = 1-minimum([1/(2*deg+1)*(v[i+1]-v[i])/((t_d[2]-t_d[1])*a), 1.0])
            else
                eta = 1-minimum([(v[i+1]-v[i])/h*1/c, 1.0])
            end
            if fluxtype == "altlr"
            # auf rho wirkende Terme (Stabilisert auf basis von outflow Term pink)
                ex_B_J0[2*i, 2*i + vshift] = -1 * eta
                ex_B_J0[2*i + 1, 2*i + vshift] = 1 * eta
                ex_B_J0E1[2*i, 2*i + vshift] = 1 * eta
                ex_B_J0E2[2*i + 1, 2*i + 1 + vshift] = -1 * eta
            # auf j wirkende Terme
                k = cellnumber + i
                im_B_J0[2*k - 1, 2*k - 1 - vshift] = epsq2 *2 * eta
                im_B_J0[2*k - 2, 2*k - 1 - vshift] = -epsq2 *2 * eta
                im_B_J0E1[2*k - 1, 2*k - 1 - vshift] = -epsq2 *2 * eta
                im_B_J0E2[2*k - 2, 2*k - 2 - vshift] = epsq2 *2 * eta
            elseif fluxtype == "full"
            # auf rho wirkende Terme
                #grün
                #ex_B[2*i-1, 2*i-1] = eps2
                ex_B_J0[2*i - 1, 2*i - 1] = -eps2 * eta
                ex_B_J0[2*i - 2, 2*i - 1] = eps2 * eta
                ex_B_J0E1[2*i - 1, 2*i - 1] = eps2 * eta
                ex_B_J0E2[2*i - 2, 2*i - 2] = -eps2 * eta
                #gelb
                #ex_B[2*i-1, 2*i-1+vshift] = -1/2
                ex_B_J0[2*i - 1, 2*i - 1 + vshift] = 1/2 * eta
                ex_B_J0[2*i - 2, 2*i - 1 + vshift] = -1/2 * eta
                ex_B_J0E1[2*i - 1, 2*i - 1 + vshift] = -1/2 * eta
                ex_B_J0E2[2*i - 2, 2*i - 2 + vshift] = 1/2 * eta                
                #pink
                #ex_B[2*i, 2*i+vshift] = 1/2
                ex_B_J0[2*i, 2*i + vshift] = -1/2 * eta
                ex_B_J0[2*i + 1, 2*i + vshift] = 1/2 * eta
                ex_B_J0E1[2*i, 2*i + vshift] = 1/2 * eta
                ex_B_J0E2[2*i + 1, 2*i + 1 + vshift] = -1/2 * eta
                # rot
                #ex_B[2*i, 2*i] = eps2
                ex_B_J0[2*i, 2*i] = -eps2 * eta
                ex_B_J0[2*i + 1, 2*i] = eps2 * eta
                ex_B_J0E1[2*i, 2*i] = eps2 * eta
                ex_B_J0E2[2*i + 1, 2*i + 1] = -eps2 * eta
            # auf j wirkende Terme
            k = cellnumber + i
                #grün
                #im_B[2*i-1, 2*i-1-vshift] = -epsq2
                im_B_J0[2*k - 1, 2*k - 1-vshift] = epsq2 * eta
                im_B_J0[2*k - 2, 2*k - 1-vshift] = -epsq2 * eta
                im_B_J0E1[2*k - 1, 2*k - 1-vshift] = -epsq2 * eta
                im_B_J0E2[2*k - 2, 2*k - 2-vshift] = epsq2 * eta
                #gelb
                #im_B[2*i-1, 2*i-1] = eps2
                im_B_J0[2*k - 1, 2*k - 1] = -eps2 * eta
                im_B_J0[2*k - 2, 2*k - 1] = eps2 * eta
                im_B_J0E1[2*k - 1, 2*k - 1] = eps2 * eta
                im_B_J0E2[2*k - 2, 2*k - 2] = -eps2 * eta
                #pink
                #im_B[2*i, 2*i] = eps2
                #ex_B[2*i, 2*i+vshift] = 1/2
                im_B_J0[2*k, 2*k] = -eps2 * eta
                im_B_J0[2*k + 1, 2*k] = eps2 * eta
                im_B_J0E1[2*k, 2*k] = eps2 * eta
                im_B_J0E2[2*k + 1, 2*k + 1] = -eps2 * eta
                #rot
                #im_B[2*i, 2*i-vshift] = epsq2
                im_B_J0[2*k, 2*k-vshift] = -epsq2 * eta
                im_B_J0[2*k + 1, 2*k-vshift] = epsq2 * eta
                im_B_J0E1[2*k, 2*k-vshift] = epsq2 * eta
                im_B_J0E2[2*k + 1, 2*k + 1-vshift] = -epsq2 * eta
            else
                throw(ArgumentError("DoD not implemented for that flux yet"))
            end
        end
    end
    return ex_B_J0, ex_B_J0E1, ex_B_J0E2, im_B_J0, im_B_J0E1, im_B_J0E2
end



function DGsemidiscretization_DoD_heat(setup, deg, nodes, num_flux; T = Float64, do_stabilize = false, c = 1.0, fix_eta = false, ext_test_func = true)
    cellnumber = setup["cellnumber"]
    v = setup["v"]
    a = setup["a"]
    u0_eval = setup["u0_eval"]
    t_d = setup["t_d"]
    h = setup["h"]
    setup["nodes"] = nodes
    setup["deg"] = deg
    setup["num_flux"] = num_flux

    splitRHS = Dict()

    # basis defines the method, basis1 will be manipulated to serve interpolation matrices for other nodes
    basis = nodes(deg, T)
    basis1 = nodes(deg, T)
    basis1_second = nodes(deg, T)

    # Some Parameters
    xrange = cellnumber*(deg+1)
    invscale(i) = sqrt(a)*2/(v[i+1]-v[i])
    shift(i) = (v[i+1]+v[i])/2
    scale(i) = (v[i+1]-v[i])/2
    
    # x_d, RHS calculation
    x_d = zeros(T, xrange)
    M = zeros(T, xrange, xrange)
    M_global = zeros(T, xrange, xrange) # TODO
    DM = zeros(T, xrange, xrange)
    DM_J1 = zeros(T, xrange, xrange)
    DM_J1_second = zeros(T, xrange, xrange)
    B = zeros(T, cellnumber*2, cellnumber*2)
    B_second = zeros(T, cellnumber*2, cellnumber*2)
    B_J0 = zeros(T, cellnumber*2, cellnumber*2)
    B_J0_second = zeros(T, cellnumber*2, cellnumber*2)
    B_J0E1 = zeros(T, cellnumber*2, cellnumber*2)
    B_J0E1_second = zeros(T, cellnumber*2, cellnumber*2)
    B_J0E2 = zeros(T, cellnumber*2, cellnumber*2)
    B_J0E2_second = zeros(T, cellnumber*2, cellnumber*2)
    V = zeros(T, cellnumber*2, xrange)
    V_J0E1 = zeros(T, cellnumber*2, xrange)
    V_J0E1_second = zeros(T, cellnumber*2, xrange)
    V_J0E2 = zeros(T, cellnumber*2, xrange)
    V_J0E2_second = zeros(T, cellnumber*2, xrange)
    TranspV = zeros(T, xrange, cellnumber*2)
    Sc = zeros(T, xrange, xrange)

    DM_J1_just_volume_cm1 = zeros(T, xrange, xrange)
    
    # Creating Transformation Matrix from reference element to real elements
    for i in 1:cellnumber
        Sc[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] = diagm(ones(T, deg+1)*invscale(i))
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
        # Constructing the Mass- and Differentiation Matrix
        M[(deg+1)*(i-1)+1:(deg+1)*i, (deg+1)*(i-1)+1:(deg+1)*i] = diagm(basis.weights)
        M_global[(deg+1)*(i-1)+1:(deg+1)*i, (deg+1)*(i-1)+1:(deg+1)*i] = diagm(basis.weights) * scale(i)
        DM[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] = basis.D' * diagm(basis.weights) 

        # Constructing the Interpolation to the vertices for the flux terms with Gauss-Legendre nodes
        #V[2*i-1, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i], basis1)
        #V[2*i, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+1], basis1)
        V[2*i-1, (deg+1)*(i-1) + 1:(deg+1)*i] = left_interpolation
        V[2*i, (deg+1)*(i-1) + 1:(deg+1)*i] = right_interpolation
        
        # Extrapolation to the Inflow of Cell + 1(V2) and the Outflow of Cut-Cell(V3) from Cell - 1 (from the perspective of Cell -1).
        if do_stabilize
            if i < cellnumber
                if 2*(v[i+2]-v[i+1])/h<1
                    V_J0E1[2*i+2, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+2], basis1) 
                    V_J0E2[2*i+3, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i+2], basis1)
                end
            end
        end

        # Extrapolation to the Inflow of Cell - 1(V2) and the Outflow of Cut-Cell(V3) from Cell + 1 (from the perspective of Cell +1).
        if do_stabilize
            if i > 1
                if 2*(v[i]-v[i-1])/h<1
                    V_J0E1_second[2*i-3, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i-1], basis1) 
                    V_J0E2_second[2*i-4, (deg+1)*(i-1) + 1:(deg+1)*i] = interpolation_matrix(v[i-1], basis1)
                end
            end
        end
            
        TranspV[(deg+1)*(i-1)+1:(deg+1)*i, 2*i-1:2*i] = V[2*i-1:2*i, (deg+1)*(i-1) + 1:(deg+1)*i]'
    end

    # Constructing numerical fluxes(with periodic BC)
    if num_flux == "Upwind" 
            for i in 1:cellnumber
                if i != 1
                    B[2*i-1, 2*(i-1)] = -1
                end
                    B[2*i, 2*i] = 1
            end
        # periodic boundary
        if setup["bcs"] == "periodic"
            B[1, end] = -1
        end
        # Upwind second term
        for i in 1:cellnumber
            B_second[2*i-1, 2*i-1] += -1
            if i != cellnumber
                B_second[2*i, 2*i+1] += 1
            end
        end
        # periodic boundary
        if setup["bcs"] == "periodic"
            B_second[end, 1] += 1
        end
    end
    

    # DoD-Stabilization
    if do_stabilize
        for i in 1:cellnumber
            if 2*(v[i+1]-v[i])/h<1
                if fix_eta == false
                    eta = 1-minimum([1/(2*deg+1)*(v[i+1]-v[i])/((t_d[2]-t_d[1])*a), 1.0])
                else
                    eta = 1-minimum([(v[i+1]-v[i])/h*1/c, 1.0])
                end
                # Reduced Outflow  from Cut-Cell (vice versa for second component)
                B_J0[2*i, 2*i] += -1 * eta
                B_J0_second[2*i - 1, 2*i - 1] += 1* eta
                # Reduced Inflow into Cell + 1 from Cut-Cell (vice versa for second component)
                B_J0[2*(i) + 1, 2*(i)] += 1 * eta
                B_J0_second[2*i - 2, 2*i-1] += -1 * eta
                # Extended Outflowflow Cut-Cell from Cell -1 (corresponds to Reduced Inflow into Cut-Cell) (vice versa for second component)
                B_J0E1[2*i , 2*i] += 1 * eta
                B_J0E1_second[2*i - 1, 2*i - 1] += -1 * eta
                # Extended Inflow into Cell + 1 from Cell - 1 (vice versa for second component)
                B_J0E2[2*(i) + 1, 2*(i) + 1]+= -1 * eta
                B_J0E2_second[2*(i) - 2, 2*(i) - 2]+= 1 * eta
                
                # Volume Term
                if deg>0
                    DM_J1[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] +=
                        - basis.D' * diagm(basis.weights) * eta
                    DM_J1_second[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] +=
                        - basis.D' * diagm(basis.weights) * eta
                    
                    basis2=nodes(deg)
                    basis2_second=nodes(deg)
                    for j in 1:deg+1
                        # Inserting these nodes into the structure basis2, to adapt the field basis1.interp_matrix to the used nodes.
                        basis2.nodes[j] = x_d[(deg+1)*(i-2) + j] 
                        basis2_second.nodes[j] = x_d[(deg+1)*(i) + j] 
                    end
                    Vol_interpol = interpolation_matrix(x_d[(deg+1)*(i-1) + 1:(deg+1)*(i)] , basis2)
                    Vol_interpol_second = interpolation_matrix(x_d[(deg+1)*(i-1) + 1:(deg+1)*(i)] , basis2_second)
                    DM_J1[(deg+1)*(i-1) + 1:(deg+1)*(i), (deg+1)*(i-2) + 1:(deg+1)*(i-1)] += 
                        + basis.D' * diagm(basis.weights) * Vol_interpol * eta
                    DM_J1_second[(deg+1)*(i-1) + 1:(deg+1)*(i), (deg+1)*(i) + 1:(deg+1)*(i+1)] += 
                        + basis.D' * diagm(basis.weights) * Vol_interpol_second * eta
                    # Extra-Stabilityterm with the extended test-function (L_E_in(u_h)-u_h)*L_E_in(v_h)
                    if ext_test_func == true
                        DM_J1[(deg+1)*(i-2) + 1:(deg+1)*(i-1), (deg+1)*(i-1) + 1:(deg+1)*i] +=
                            Vol_interpol' *  basis.D' * diagm(basis.weights) * eta
                        DM_J1_second[(deg+1)*(i) + 1:(deg+1)*(i+1), (deg+1)*(i-1) + 1:(deg+1)*i] +=
                            Vol_interpol_second' *  basis.D' * diagm(basis.weights) * eta
                        DM_J1[(deg+1)*(i-2) + 1:(deg+1)*(i-1), (deg+1)*(i-2) + 1:(deg+1)*(i-1)] += 
                            - Vol_interpol' * basis.D' * diagm(basis.weights) * Vol_interpol * eta
                        DM_J1_second[(deg+1)*(i) + 1:(deg+1)*(i+1), (deg+1)*(i) + 1:(deg+1)*(i+1)] += 
                            - Vol_interpol_second' * basis.D' * diagm(basis.weights) * Vol_interpol_second * eta
                        DM_J1_just_volume_cm1[(deg+1)*(i-2) + 1:(deg+1)*(i-1), (deg+1)*(i-2) + 1:(deg+1)*(i-1)] = DM_J1[(deg+1)*(i-2) + 1:(deg+1)*(i-1), (deg+1)*(i-2) + 1:(deg+1)*(i-1)]
                    end
                end
            end
        end
    end
    
    # Construncting RHS Matrix
    invM = inv(M)
    # DoD: linear advection:
    # RHS_mat = Sc*invM * (-TranspV* (B*V + B_J0*V + B_J0E1*V_J0E1 + B_J0E2*V_J0E2) + DM + DM_J1)
    # heat without DoD
    RHS_mat = Sc*invM * ((-TranspV* B*V + DM)*Sc*invM*(-TranspV* B_second*V + DM))
    # heat with DoD
    RHS_mat = Sc*invM * ((-TranspV* (B*V + B_J0*V + B_J0E1*V_J0E1 + B_J0E2*V_J0E2) + DM + DM_J1)*Sc*invM*(-TranspV* (B_second*V + B_J0_second*V + B_J0E1_second*V_J0E1_second + B_J0E2_second*V_J0E2_second) + DM + DM_J1_second))
    
    #display(Sc*invM * ((-TranspV* (B*V + B_J0*V + B_J0E1*V_J0E1 + B_J0E2*V_J0E2) + DM + DM_J1)))

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
    splitRHS["DM_J1_just_vol_cm1"] = Sc*invM * DM_J1_just_volume_cm1
    # Constructing the initial conition
    u0 = [u0_eval(x) for x in x_d]
    setup["u0"] = u0
    setup["M"] = M_global
    setup["x_d"] = x_d
    setup["basis"] = basis
    
    #display((-TranspV* (B*V + B2*V2))[1:9, 1:9])
    #display(DM_J1[1:10, 1:10])
    
    return RHS_mat, setup, splitRHS
end
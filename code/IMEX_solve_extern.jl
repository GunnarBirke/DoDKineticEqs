function IMEX_solve_extern(problem; TMM = IMEXEuler)

    ex_RHS_mat = problem["ex_RHS_mat"]
    im_RHS_mat = problem["im_RHS_mat"]
    Tmax = problem["Tmax"]
    tspan = [0.0,Tmax];
    x_d = problem["x_d"];
    u0 = problem["u0"]
    N = problem["cellnumber"]
    CFL = problem["CFL"]

    function A_ex!(du, u, p, t)
        du .= ex_RHS_mat*u
    end
    function A_im!(du, u, p, t)
        du .= im_RHS_mat*u
    end
    
    
    f = SplitFunction(A_im!, A_ex!)
    prob = SplitODEProblem(f,u0, tspan)
    dt = 2*pi*CFL/N^2
    sol = solve(prob, TMM(), dt = dt)
    problem["sol"] = sol
    (xsteps, tsteps) = size(sol)
    
    # Telegraph exact solution (with u_0_rho = 1/rsinx, u_0_j=cos(x))
    r = -2/(1+sqrt(1-4*epsilon^2))
    u_exact = zeros((deg+1)*2*N, tsteps);
    u_exact[1:(deg+1)*N, 1:tsteps] = [1/r*exp(r*t)*sin(x) for x in x_d, t in range(0.0, Tmax, length = tsteps)];
    u_exact[(deg+1)*N+1:(deg+1)*2*N, 1:tsteps] = [exp(r*t)*cos(x) for x in x_d, t in range(0.0, Tmax, length = tsteps)];
    problem["u_exact"] = u_exact

    v = problem["v"]
    nodes = problem["nodes"]
    sqerror = 0
    for i in 1:problem["cellnumber"]
        basis1 = nodes(deg)
        for j in 1:deg+1
            # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
            basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
        end
        sqerror += integrate((sol[(deg+1)*(i-1)+1:(deg+1)*i, end] - u_exact[(deg+1)*(i-1)+1:(deg+1)*i, end]).^2, basis1.weights.*(v[i+1]-v[i])/2)
    end
    problem["error"] = sqrt(sqerror)

    return problem
end
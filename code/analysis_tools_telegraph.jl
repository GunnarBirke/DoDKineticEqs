#############################################################################
############################ Convergence test ###############################
#############################################################################

function test_convergence_telegraph_background(basis, deg ; Tmax = 5.0, TMM = SSPRK3, exprange = [3,8], CFL = 0.1, a=1, epsilon = 0.4)
    stepsizes = 2 .^range(exprange[1],exprange[2])
    errors = zeros(length(stepsizes))
    for (istep, stepsize) in enumerate(stepsizes)
        problem = setup_problem_eq("telsin", -pi, pi, stepsize, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = epsilon);
        RHS_mat, problem = DGsemidiscretization_DoD_telegraph(problem, deg, basis, "Upwind",  do_stabilize = false, eq_type = "telegraph", fluxtype = "full");
        solution, u_exact = TMM(problem, RHS_mat)
        #errors[istep] = norm(solution[:, end] - u_exact[:, end])
        #println(solution[:, end])

        x_d = problem["x_d"]
        v = problem["v"]
        deg = problem["deg"]
        nodes = problem["nodes"]
        errors[istep] = 0
        for i in 1:problem["cellnumber"]
            basis1 = nodes(deg)
            for j in 1:deg+1
                # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
                basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
            end
            errors[istep] += integrate((solution[(deg+1)*(i-1)+1:(deg+1)*i, end] - u_exact[(deg+1)*(i-1)+1:(deg+1)*i, end]).^2, basis1.weights.*(v[i+1]-v[i])/2)
        end
        errors[istep] = sqrt(errors[istep])
        
    end
    return stepsizes, errors
end
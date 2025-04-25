using LinearAlgebra, Plots, PolynomialBases, OrdinaryDiffEq, PrettyTables, DelimitedFiles
include("../functions.jl")
include("../vizualize.jl")
include("../functionsRK.jl")
include("../IMEX_solve_extern.jl")


function upwindSBP_test( ;
                            N = 2^4, 
                            deg = 0, 
                            epsilon = 0.5, 
                            basis = GaussLegendre, 
                            c = 0.5, 
                            ext_test_func = true,
                            digit_tolerance = 12, 
                            alphas = [0.4, 0.3, 0.1, 10^-3, 10^-8],
                            #### just relevant for plot
                            plot_sol = false,
                            TMM = IMEXEuler, # Just for IMEX extern solver needed
                            CFL = 1/2*(deg+1)*epsilon/N,
                            Tmax = 3.0,
                            print_result = true,
                            )

    problem = setup_problem_eq("telsin", -pi, pi, N, Tmax = Tmax, a = 1.0, CFL = CFL, bcs = "periodic", epsilon = epsilon);
    for (ialpha, alpha) in enumerate(alphas)
        problem = include_cut_cell(problem, alpha, 3 + 3*(ialpha-1));
    end
    #### discretize in space
    RHS_mat, problem, splitRHS = DGsemidiscretization_DoD_telegraph(problem, deg, basis, "Upwind",  do_stabilize = true, fix_eta = true, c = c, fluxtype = "altlr", ext_test_func = ext_test_func);
    ex_RHS_mat = problem["ex_RHS_mat"]
    im_RHS_mat = problem["im_RHS_mat"]
    eigenvalues = eigvals(splitRHS["symm_negsemi"])
    generalized_eigenvalues = eigvals(1/2*(splitRHS["symm_negsemi"] + splitRHS["symm_negsemi"]'))

    # Upwind-Operators
    if print_result == true
        isupwindsbp = true
        println("Results:")
        if round.(splitRHS["IntByParts"], digits = digit_tolerance ) != zeros(size(splitRHS["IntByParts"]))
            isupwindspb = false
            println("-IntByParts not fulfilled")
        end
        if round.(splitRHS["symm_negsemi"], digits = digit_tolerance ) != round.(splitRHS["symm_negsemi"]', digits = digit_tolerance )
            isupwindsbp = false
            println("-Necessary operator not symmetric")
        end
        for eigval in eigenvalues
            if round(imag(eigval), digits = digit_tolerance ) != 0 || round.(real(eigval), digits = digit_tolerance )>0
                isupwindsbp = false
                println("-Wrong eigenvalue: $(eigval)")
            end
        end
        if isupwindsbp
            println("-> SBP Upwind confirmed!")
        else
            println("-> SBP Upwind disproven!")
        end
        println("------------------------------------------------------------------------------------")
        isgenupwindsbp = true
        if round.(splitRHS["IntByParts"], digits = digit_tolerance ) != -round.(splitRHS["IntByParts"], digits = digit_tolerance )'
            isgenupwindsbp = false
            println("-generalized IntByParts not fulfilled")
        end
        if round.(1/2*(splitRHS["symm_negsemi"] + splitRHS["symm_negsemi"]'), digits = digit_tolerance ) != round.(1/2*(splitRHS["symm_negsemi"] + splitRHS["symm_negsemi"]')', digits = digit_tolerance )
            isgenupwindsbp = false
            println("-Necessary symmetric operator not really symmetric")
        end
        for eigval in generalized_eigenvalues
            if round(imag(eigval), digits = digit_tolerance ) != 0 || round.(real(eigval), digits = digit_tolerance )>0
                isgenupwindsbp = false
                println("-Wrong eigenvalue of symmetric operator: $(eigval)")
            end
        end
        if isgenupwindsbp
            println("-> generalized SBP Upwind confirmed!")
        else
            println("-> generalized SBP Upwind disproven!")
        end
    end

    if plot_sol == true
        #sol, u_exact = SSPRK3(problem, RHS_mat)
        sol, u_exact = ImExEuler(problem, RHS_mat)
        #sol, u_exact = ARS2(problem, RHS_mat)
        #sol, u_exact = ARS3(problem, RHS_mat)

        full_sol = sol[1:Int(size(sol)[1]/2), :] .+sol[Int(size(sol)[1]/2+1):size(sol)[1], :]*epsilon
        
        t_dest = size(sol)[2] # Observing endpoint
        x_d = problem["x_d"];
        plot(x_d, u_exact[1:Int(size(u_exact)[1]/2), t_dest] .+ epsilon*u_exact[Int(size(u_exact)[1]/2+1):size(u_exact)[1], t_dest], label = "exact")
        plot!(x_d, full_sol[:, t_dest], linestyle = :dash, label = "approx")
    end
    return splitRHS["IntByParts"], splitRHS["symm_negsemi"], eigenvalues, problem, splitRHS
end

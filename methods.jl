using JuMP
using Gurobi

function convex_hull_approx()

    ## ----------- DATA SAMPLING FOR CONVEX HULL APPROXIMATION ----------- ##

    # Independent Sample Vectors
    V1_sample = range(min_V1, max_V1, length = M)
    u1_sample = range(min_ut1, max_ut1, length = N)

    V2_sample = range(min_V2, max_V2, length = M)
    u2_sample = range(min_ut2, max_ut2, length = N)

    # Define Hydraulic Head Function
    f1(V) = a1 * (V^b1)
    f2(V) = a2 * (V^b2)

    # Precomputed Matrix of Power Outputs
    P1 = [eta * g * rho_w * f1(V) * u for u in u1_sample, V in V1_sample] / (3.6e9)
    P2 = [eta * g * rho_w * f2(V) * u for u in u2_sample, V in V2_sample] / (3.6e9)
    
    # Create the optimization model
    model = Model(Gurobi.Optimizer)

    ## ----------- SIMULATIONS ----------- ##

    ### Unit 01: Bonneville Dam (Downstream)
    ## Define variables
    @variable(model, s1[1:T] >= 0) # Spill Outflow [m3/hr]
    @variable(model, lam1[1:N, 1:M, 1:T] >= 0)  # lambda matrix

    ## Initial conditions
    @constraint(model, MassBalInit1, sum(lam1[i,j,1] * V1_sample[j] for i=1:N, j=1:M) == V0_01)
    @constraint(model, RampRateInit1, sum(lam1[i,j,1] * u1_sample[i] for i=1:N, j=1:M) == min_ut1)

    ## Expressions
    @expression(model, V1[t=1:T], sum(lam1[i,j,t] * V1_sample[j] for i=1:N, j=1:M))
    @expression(model, u1[t=1:T], sum(lam1[i,j,t] * u1_sample[i] for i=1:N, j=1:M))
    @expression(model, p1[t=1:T], sum(lam1[i,j,t] * P1[i,j]      for i=1:N, j=1:M))

    ## Constraints
    @constraint(model, MassBal1[t in 2:T], V1[t] == V1[t-1] + q1[t] - u1[t] - s1[t])
    @constraint(model, Release1[t in 2:T], min_ut1 <= u1[t] <= max_ut1)
    @constraint(model, RampRate1[t in 2:T], RR_dn1 <= u1[t] - u1[t-1] <= RR_up1)
    @constraint(model, FeederCap1[t in 1:T], 0 <= p1[t] <= F1)
    @constraint(model, Volume1[t in 1:T], min_V1 <= V1[t] <= max_V1)
    @constraint(model, lambda1[t in 1:T], sum(lam1[:,:,t]) == 1)

    ### Unit 02: Dalles Dam (Downstream)
    ## Define variables
    @variable(model, s2[1:T] >= 0) # Spill Outflow [m3/hr]
    @variable(model, lam2[1:N, 1:M, 1:T] >= 0)  # lambda matrix

    ## Initial conditions
    @constraint(model, MassBalInit2, sum(lam2[i,j,1] * V2_sample[j] for i=1:N, j=1:M) == V0_02)
    @constraint(model, RampRateInit2, sum(lam2[i,j,1] * u2_sample[i] for i=1:N, j=1:M) == min_ut2)

    ## Expressions
    @expression(model, V2[t=1:T], sum(lam2[i,j,t] * V2_sample[j] for i=1:N, j=1:M))
    @expression(model, u2[t=1:T], sum(lam2[i,j,t] * u2_sample[i] for i=1:N, j=1:M))
    @expression(model, p2[t=1:T], sum(lam2[i,j,t] * P2[i,j]      for i=1:N, j=1:M))

    ## Constraints
    @constraint(model, MassBal2[t in 2:T], V2[t] == V2[t-1] + q2[t] - u2[t] - s2[t])
    @constraint(model, Release2[t in 2:T], min_ut2 <= u2[t] <= max_ut2)
    @constraint(model, RampRate2[t in 2:T], RR_dn2 <= u2[t] - u2[t-1] <= RR_up2)
    @constraint(model, FeederCap2[t in 1:T], 0 <= p2[t] <= F2)
    @constraint(model, Volume2[t in 1:T], min_V2 <= V2[t] <= max_V2)
    @constraint(model, lambda2[t in 1:T], sum(lam2[:,:,t]) == 1)

    # Objective function
    @objective(model, Max, sum(p1 + p2)) 

    # Solve the optimization problem
    optimize!(model)

    # Total Generation
    obj = objective_value(model);

    return model, obj, value.(s1), value.(lam1), value.(V1), value.(u1), value.(p1), value.(s2), value.(lam2), value.(V2), value.(u2), value.(p2)

end


function MINLP()
    
    # Create the optimization model
    model = Model(Gurobi.Optimizer)

    ## Define variables
    # Unit 01: Bonneville Dam (Downstream)
    @variable(model, V1[1:T] >= 0) # Reservoir Volume [m3]
    @variable(model, p1[1:T] >= 0) # Generation [MWh]
    @variable(model, u1[1:T])      # Generation Outflow [m3/hr]
    @variable(model, s1[1:T] >= 0) # Spill Outflow [m3/hr]

    # Unit 02: The Dalles Dam (Upstream)
    @variable(model, V2[1:T] >= 0) # Reservoir Volume [m3]
    @variable(model, p2[1:T] >= 0) # Generation [MWh]
    @variable(model, u2[1:T])      # Generation Outflow [m3/hr]
    @variable(model, s2[1:T] >= 0) # Spill Outflow [m3/hr]

    ## Initial conditions
    # Unit 1
    @constraint(model, MassBalInit1, V1[1] == V0_01)
    @constraint(model, RampRateInit1, u1[1] == min_ut1)
    # Unit 2
    @constraint(model, MassBalInit2, V2[1] == V0_02)
    @constraint(model, RampRateInit2, u2[1] == min_ut2)

    # Objective function
    @objective(model, Max, sum(p1 + p2)) 

    ## Constraints
    # Unit 1
    @constraint(model, MassBal1[t in 2:T], V1[t] == V1[t-1] + q1[t] - u1[t] - s1[t])
    @constraint(model, ReleaseEnergy1[t in 1:T], p1[t] <= (eta * g * rho_w * u1[t] * a1 * (V1[t]^b1))/(3.6e9))
    @constraint(model, Release1[t in 2:T], min_ut1 <= u1[t] <= max_ut1)
    @constraint(model, RampRate1[t in 2:T], RR_dn1 <= u1[t] - u1[t-1] <= RR_up1)
    @constraint(model, FeederCap1[t in 1:T], 0 <= p1[t] <= F1)
    @constraint(model, Volume1[t in 1:T], min_V1 <= V1[t] <= max_V1)

    # Unit 2
    @constraint(model, MassBal2[t in 2:T], V2[t] == V2[t-1] + q2[t] - u2[t] - s2[t])
    @constraint(model, ReleaseEnergy2[t in 1:T], p2[t] <= (eta * g * rho_w * u2[t] * a2 * (V2[t]^b2))/(3.6e9))
    @constraint(model, Release2[t in 2:T], min_ut2 <= u2[t] <= max_ut2)
    @constraint(model, RampRate2[t in 2:T], RR_dn2 <= u2[t] - u2[t-1] <= RR_up2)
    @constraint(model, FeederCap2[t in 1:T], 0 <= p2[t] <= F2)
    @constraint(model, Volume2[t in 1:T], min_V2 <= V2[t] <= max_V2)

    #= Set parameters for piecewise-linear approximation (FASTER)
    set_attribute(model, "FuncNonlinear", 0)    # Use piecewise-linear (default)
    set_attribute(model, "FuncPieces", 4)       
    set_attribute(model, "FuncPieceError", 5e-2)  # Allow larger errors
    =#   

    # Solve the optimization problem
    optimize!(model)

    # Revenue
    obj = objective_value(model)

    return model, obj, value.(s1), value.(V1), value.(u1), value.(p1), value.(s2), value.(V2), value.(u2), value.(p2)

end


function model_report(model)
    # Check the solution method used
    method_used = MOI.get(model, Gurobi.ModelAttribute("ConcurrentWinMethod"))
    status = MOI.get(model, MOI.TerminationStatus())

    # Check iteration counts
    iter_count = MOI.get(model, Gurobi.ModelAttribute("IterCount"))
    bar_iter_count = MOI.get(model, Gurobi.ModelAttribute("BarIterCount"))

    # Check if there are any general constraints (which could make it behave like MIP)
    num_gen_constrs = MOI.get(model, Gurobi.ModelAttribute("NumGenConstrs"))
    num_sos = MOI.get(model, Gurobi.ModelAttribute("NumSOS"))

    println()
    println("--- MODEL REPORT ---")
    println("Method used: ", method_used)  # 0=primal simplex, 1=dual simplex, 2=barrier
    println("Status: ", status)
    println("Simplex iterations: ", iter_count)
    println("Barrier iterations: ", bar_iter_count)
    println("General constraints: ", num_gen_constrs)
    println("SOS constraints: ", num_sos)
end

function variable_report(method, obj, p1, u1, s1, p2, u2, s2)

    println()
    println("--- VARIABLE REPORT ---")
    println("Method used: " * method)  
    if method == "CHA"
        println("M = ", M)
        println("N = ", N)
    end 
    println("Total Generation (obj): ", round(obj, sigdigits=3))
    println("---Cumulative Values---")
    println("Generation 01 [MWh]: ", round(sum(p1), sigdigits=3))
    println("Generation 02 [MWh]: ", round(sum(p2), sigdigits=3))
    println("Generation Release 01 [m3]: ", round(sum(u1), sigdigits=3))
    println("Generation Release 02 [m3]: ", round(sum(u2), sigdigits=3))
    println("Spill Release 01 [m3]: ", round(sum(s1), sigdigits=3))
    println("Spill Release 02 [m3]: ", round(sum(s2), sigdigits=3))
    println()

end
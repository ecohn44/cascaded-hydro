using JuMP
using Gurobi
using Ipopt


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


function pred_inflow_diu()

end


function MINLP()
    
    # Create the optimization model
    # model = Model(Gurobi.Optimizer)
    model = Model(Ipopt.Optimizer)

    ## Define variables
    # Unit 01: Bonneville Dam (Downstream)
    @variable(model, V1[1:T] >= 0) # Reservoir Volume [m3]
    @variable(model, p1[1:T] >= 0) # Generation [MWh]
    @variable(model, u1[1:T] >= 0) # Generation Outflow [m3/hr]
    @variable(model, s1[1:T] >= 0) # Spill Outflow [m3/hr]

    # Unit 02: The Dalles Dam (Upstream)
    @variable(model, V2[1:T] >= 0) # Reservoir Volume [m3]
    @variable(model, p2[1:T] >= 0) # Generation [MWh]
    @variable(model, u2[1:T] >= 0) # Generation Outflow [m3/hr]
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
    @constraint(model, ReleaseEnergy1[t in 1:T], p1[t] == (eta * g * rho_w * u1[t] * a1 * (V1[t]^b1))/(3.6e9))
    @constraint(model, Release1[t in 2:T], min_ut1 <= u1[t] <= max_ut1)
    @constraint(model, RampRate1[t in 2:T], RR_dn1 <= u1[t] - u1[t-1] <= RR_up1)
    @constraint(model, FeederCap1[t in 1:T], 0 <= p1[t] <= F1)
    @constraint(model, Volume1[t in 1:T], min_V1 <= V1[t] <= max_V1)

    # Unit 2
    @constraint(model, MassBal2[t in 2:T], V2[t] == V2[t-1] + q2[t] - u2[t] - s2[t])
    @constraint(model, ReleaseEnergy2[t in 1:T], p2[t] == (eta * g * rho_w * u2[t] * a2 * (V2[t]^b2))/(3.6e9))
    @constraint(model, Release2[t in 2:T], min_ut2 <= u2[t] <= max_ut2)
    @constraint(model, RampRate2[t in 2:T], RR_dn2 <= u2[t] - u2[t-1] <= RR_up2)
    @constraint(model, FeederCap2[t in 1:T], 0 <= p2[t] <= F2)
    @constraint(model, Volume2[t in 1:T], min_V2 <= V2[t] <= max_V2)

    # Solve the optimization problem
    optimize!(model)

    # Revenue
    obj = objective_value(model)

    return model, obj, value.(s1), value.(V1), value.(u1), value.(p1), value.(s2), value.(V2), value.(u2), value.(p2)

end


function MINLP_loop(q1_s, q2_s) #, framework)
    
    # Initialize empty vectors to store DVs
    V1_s = zeros(Float64, T)
    p1_s = zeros(Float64, T)
    u1_s = zeros(Float64, T)
    s1_s = zeros(Float64, T)
    q1_pred = zeros(Float64, T)

    V2_s = zeros(Float64, T)
    p2_s = zeros(Float64, T)
    u2_s = zeros(Float64, T)
    s2_s = zeros(Float64, T)

    # Create the optimization model
    model = Model(Gurobi.Optimizer)
    # model = Model(Ipopt.Optimizer)

    ## Define variables
    # Unit 01: Bonneville Dam (Downstream)
    @variable(model, V1 >= 0) # Reservoir Volume [m3]
    @variable(model, p1 >= 0) # Generation [MWh]
    @variable(model, u1 >= 0) # Generation Outflow [m3/hr]
    @variable(model, s1 >= 0) # Spill Outflow [m3/hr]

    # Unit 02: The Dalles Dam (Upstream)
    @variable(model, V2 >= 0) # Reservoir Volume [m3]
    @variable(model, p2 >= 0) # Generation [MWh]
    @variable(model, u2 >= 0) # Generation Outflow [m3/hr]
    @variable(model, s2 >= 0) # Spill Outflow [m3/hr]

    # Objective function
    @objective(model, Max, sum(p1 + p2))

    ## Constraints
    # Unit 1
    @constraint(model, MassBal1, V1 + u1 + s1 == V1_s[1] + q1_s[1])
    @constraint(model, ReleaseEnergy1, p1 == (eta * g * rho_w * u1 * a1 * (V1^b1))/(3.6e9))
    @constraint(model, Release1, min_ut1 <= u1 <= max_ut1)
    @constraint(model, RampRateDn1,  -u1 <= -(RR_dn1 + u1_s[1]))
    @constraint(model, RampRateUp1, u1 <= RR_up1 + u1_s[1])
    @constraint(model, FeederCap1, 0 <= p1 <= F1)
    @constraint(model, Volume1, min_V1 <= V1 <= max_V1)

    # Unit 2
    @constraint(model, MassBal2, V2 + u2 + s2 == V2_s[1] + q2_s[1])
    @constraint(model, ReleaseEnergy2, p2 == (eta * g * rho_w * u2 * a2 * (V2^b2))/(3.6e9))
    @constraint(model, Release2, min_ut2 <= u2 <= max_ut2)
    @constraint(model, RampRateDn2,  -u2 <= -(RR_dn2 + u2_s[1]))
    @constraint(model, RampRateUp2, u2 <= RR_up2 + u2_s[1])
    @constraint(model, FeederCap2, 0 <= p2 <= F2)
    @constraint(model, Volume2, min_V2 <= V2 <= max_V2) 

    for t = 1:T

        #=
        ## Predict inflow at current time step 
        if framework == "DET"
            q1_pred[t] = q1_s[t]
        elseif framework == "DIU"
            q1_pred[t] = q1_s[t]
        end
        =#
        q1_pred[t] = q1_s[t]

        if t == 1
            ## Set initial flow conditions 
            # Unit 1
            set_normalized_rhs(RampRateUp1, min_ut1)
            set_normalized_rhs(RampRateDn1, -min_ut1)
            # Unit 2
            set_normalized_rhs(RampRateUp2, min_ut2)
            set_normalized_rhs(RampRateDn2, -min_ut2)
            
            ## Set intial volume conditions
            set_normalized_rhs(MassBal1, V0_01 + q1_pred[1])
            set_normalized_rhs(MassBal2, V0_02 + q2_s[1])
        else
            # Update ramp rate with previous flow
            # Unit 1
            set_normalized_rhs(RampRateUp1, RR_up1 + u1_s[t-1])
            set_normalized_rhs(RampRateDn1, -(RR_dn1 + u1_s[t-1]))
            # Unit 2
            set_normalized_rhs(RampRateUp2, RR_up2 + u2_s[t-1])
            set_normalized_rhs(RampRateDn2, -(RR_dn2 + u2_s[t-1]))

            # Update mass balance with previous volume
            set_normalized_rhs(MassBal1, V1_s[t-1] + q1_pred[t])
            set_normalized_rhs(MassBal2, V2_s[t-1] + q2_s[t])
        end
 
        # Solve the optimization problem
        optimize!(model)

        # Store decision variables
        V1_s[t] = value.(V1)
        p1_s[t] = value.(p1)
        u1_s[t] = value.(u1)
        s1_s[t] = value.(s1)
        V2_s[t] = value.(V2)
        p2_s[t] = value.(p2)
        u2_s[t] = value.(u2)
        s2_s[t] = value.(s2)

    end 

    obj = sum(p1_s + p2_s)

    return model, obj, V1_s, p1_s, u1_s, s1_s, q1_pred, V2_s, p2_s, u2_s, s2_s

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

    sd = 4

    println()
    println("--- VARIABLE REPORT ---")
    println("Method used: " * method)  
    if method == "CHA"
        println("M = ", M)
        println("N = ", N)
    end 
    println("Total Generation (obj): ", round(obj, sigdigits=sd+1))
    println("---Cumulative Values---")
    println("Generation 01 [MWh]: ", round(sum(p1), sigdigits=sd+1))
    println("Generation 02 [MWh]: ", round(sum(p2), sigdigits=sd+1))
    println("Generation Release 01 [m3]: ", round(sum(u1), sigdigits=sd))
    println("Generation Release 02 [m3]: ", round(sum(u2), sigdigits=sd))
    println("Spill Release 01 [m3]: ", round(sum(s1), sigdigits=sd-1))
    println("Spill Release 02 [m3]: ", round(sum(s2), sigdigits=sd-1))
    println()

end


function approximation_error(u, V, y, b)
    for t in 1:T
        # True nonlinear value
        true_y = value(u[t]) * (value(V[t]))^b
        
        # Approximated value
        approx_y = value(y[t])
        
        # Error
        error = abs(true_y - approx_y)

        # Percent Error
        percent_error = (error / abs(true_y)) * 100
        
        println("t=$t: True y1=$true_y, Approx y1=$approx_y, Error=$error")
        println("t=$t: Percent Error=$percent_error")
        
        if percent_error > 1e-3  # Flag significant errors
            println("  *** SIGNIFICANT APPROXIMATION ERROR ***")
        end
        println()
    end
end


function [model, obj, X, std_hat, phi_vals, alpha_vals, U_eff] = optimization(T, N, c, q, lag, framework, bounds, params, s)

    % Initialize decision variable storage
    % X columns: 1=V1, 2=p1, 3=u1, 4=s1, 5=q1, 
    %            6=V2, 7=p2, 8=u2, 9=s2, 10=q2
    X = zeros(T,10); % 5*n
    std_hat = zeros(T,2); 
    alpha_vals = zeros(T,2); 
    phi_vals = zeros(T);
    U_eff = zeros(T,4); % 2*n

    % Extract streamflow time series 
    q1_s = q(:,1); % historical reference for upstream unit (already lagged)

    % Define risk level 
    eps_t = 0.99;
    n = 2; % number of units 

    % Define decision variables (YALMIP)
    yalmip('clear');

    % Unit 01
    V1 = sdpvar(T,1);
    p1 = sdpvar(T,1);
    u1 = sdpvar(T,1);
    s1 = sdpvar(T,1);

    % Unit 02
    V2 = sdpvar(T,1);
    p2 = sdpvar(T,1);
    u2 = sdpvar(T,1);
    s2 = sdpvar(T,1);

    % Objective
    cons = [];
    Objective = sum((p1 + p2) - (s1 + s2));

    %% Static Constraints
    % Water Release Bounds
    cons = [cons, s(1).min_ut <= u1 <= s(1).max_ut];
    cons = [cons, s(2).min_ut <= u2 <= s(2).max_ut];

    % Feeder Capacity
    cons = [cons, 0 <= p1 <= s(1).F];
    cons = [cons, 0 <= p2 <= s(2).F];
    
    % Volume Bounds
    cons = [cons, s(1).min_V <= V1 <= s(1).max_V];
    cons = [cons, s(2).min_V <= V2 <= s(2).max_V];

    % Spill Bounds
    cons = [cons, 0 <= s1];
    cons = [cons, 0 <= s2];


    %% Build constraints
    for t = 1:T
       
        %% Forecast Inflow
        % q_i = predicted mean inflow, stdhat_i = prdicted std dev of inflow 
        [q1, q2, std_hat(t,:)] = forecast_inflow(X, t, q1_s, lag, framework, params, s);

        %% Determine Volume Bounds
        switch bounds
            case {"det", "jcc-ssh"}
                % Deterministic 
                q1_min = q1;
                q1_max = q1; 
                q2_min = q2;
                q2_max = q2;

            case "icc"   
                % Calculate z-score 
                z = norminv(eps_t); 

                % Individual Chance Constraints
                q1_min = q1 - z*std_hat(t,1);
                q1_max = q1 + z*std_hat(t,1);
                q2_min = q2 - z*std_hat(t,2);
                q2_max = q2 + z*std_hat(t,2);

            case "jcc-bon"
                % Calculate z-score 
                z = norminv(1 - (1 - eps_t)/n);

                % Individual Chance Constraints
                q1_min = q1 - z*std_hat(t,1);
                q1_max = q1 + z*std_hat(t,1);
                q2_min = q2 - z*std_hat(t,2);
                q2_max = q2 + z*std_hat(t,2);
        end

        % Compute allowable outflow bands (m3/hr) 
        if t == 1
            Vprev1 = s(1).V0;  Vprev2 = s(2).V0;
        else
            Vprev1 = X(t-1,1);  Vprev2 = X(t-1,6);
        end

        % Storage-implied band from chance constraints
        flow_max1 = Vprev1 + q1_min - s(1).min_V;   % ensures V >= min_V
        flow_min1 = Vprev1 + q1_max - s(1).max_V;   % ensures V <= max_V
        flow_max2 = Vprev2 + q2_min - s(2).min_V;
        flow_min2 = Vprev2 + q2_max - s(2).max_V;

        % Intersect chance constrained flow bounds with ramps
        if t == 1
            u1_lo = s(1).min_ut;  u2_lo = s(2).min_ut;
            u1_hi = s(1).min_ut;  u2_hi = s(2).min_ut;
        else
            % signed convention: RR_dn <= 0, RR_up >= 0
            u1_lo = max(s(1).min_ut, X(t-1,3) + s(1).RR_dn);
            u2_lo = max(s(2).min_ut, X(t-1,8) + s(2).RR_dn);
            u1_hi = min(s(1).max_ut, X(t-1,3) + s(1).RR_up);
            u2_hi = min(s(2).max_ut, X(t-1,8) + s(2).RR_up);
        end
        
        % Clip both sides to [u_lo, u_hi] so the band is truly feasible
        flow_max1 = min(max(flow_max1, u1_lo), u1_hi);
        flow_min1 = min(max(flow_min1, u1_lo), u1_hi);
        flow_max2 = min(max(flow_max2, u2_lo), u2_hi);
        flow_min2 = min(max(flow_min2, u2_lo), u2_hi);
        
        % Keep ordering (in case clipping collapses interval)
        if flow_min1 > flow_max1, flow_min1 = flow_max1; end
        if flow_min2 > flow_max2, flow_min2 = flow_max2; end
        
        % Store effective flow bounds
        U_eff(t,:) = [flow_max1, flow_min1, flow_max2, flow_min2];


        %% Time-Varying Constraints
        if t == 1 % Initial conditions
            % Volume Bounds
            cons = [cons, u1(t) + s1(t) <= flow_max1];   
            cons = [cons, u1(t) + s1(t) >= q1_max + s(1).V0 - s(1).max_V];
            cons = [cons, u2(t) + s2(t) <= flow_max2];   
            cons = [cons, u2(t) + s2(t) >= q2_max + s(2).V0 - s(2).max_V];  


            % Mass Balance
            cons = [cons, V1(t) + u1(t) + s1(t) == s(1).V0 + q1];
            cons = [cons, V2(t) + u2(t) + s2(t) == s(2).V0 + q2];
            
            % Ramp Rates
            cons = [cons, u1(t) == s(1).min_ut];
            cons = [cons, u2(t) == s(2).min_ut];

            % Power Production
            cons = [cons, p1(t) == c*u1(t)*s(1).a*(s(1).V0^s(1).b)];
            cons = [cons, p2(t) == c*u2(t)*s(2).a*(s(2).V0^s(2).b)];
        
        else

            % Mass Balance
            cons = [cons, V1(t) + u1(t) + s1(t) == X(t-1,1) + q1];
            cons = [cons, V2(t) + u2(t) + s2(t) == X(t-1,6) + q2];
            
            % Ramp Rate Up
            cons = [cons, u1(t) <= s(1).RR_up + u1(t-1)];
            cons = [cons, u2(t) <= s(2).RR_up + u2(t-1)];

            % Ramp Rate Down
            cons = [cons, -u1(t) <= -(s(1).RR_dn + u1(t-1))];
            cons = [cons, -u2(t) <= -(s(2).RR_dn + u2(t-1))];

            % Power Production
            h1_lookup = s(1).a * X(t-1,1)^s(1).b; % previous hydraulic head
            h1_ref = map_V_to_h(h1_lookup, s(1).h_lbounds, s(1).h_rbounds, s(1).h_refvals);
            cons = [cons, p1(t) == c*u1(t)*h1_ref];

            h2_lookup = s(2).a * X(t-1,6)^s(2).b; % previous hydraulic head
            h2_ref = map_V_to_h(h2_lookup, s(2).h_lbounds, s(2).h_rbounds, s(2).h_refvals);
            cons = [cons, p2(t) == c*u2(t)*h2_ref];

            switch bounds
                case {"det","icc","jcc-bon"}
                    % Add deterministic / ICC / Bonferroni volume bounds
                    cons = [cons, u1(t) + s1(t) <= flow_max1];  
                    cons = [cons, u1(t) + s1(t) >= X(t-1,1) + q1_max - s(1).max_V];
                    cons = [cons, u2(t) + s2(t) <= flow_max2];  
                    cons = [cons, u2(t) + s2(t) >= X(t-1,6) + q2_max - s(2).max_V]; 
            end
                          
        end

        %% Solve Optimization Problem 
        options = sdpsettings('solver','gurobi','verbose',0);
        switch bounds
            case {"det","icc","jcc-bon"}
                % Standard solve
                model = optimize(cons, -Objective, options);

                % Abort early on infeasible/numerical issues
                if model.problem ~= 0
                    error('Solver issue at t=%d: %s', t, yalmiperror(model.problem));
                end
                    
                X(t,1) = value(V1(t));
                X(t,2) = value(p1(t));
                X(t,3) = value(u1(t));
                X(t,4) = value(s1(t));
                X(t,5) = q1;
                X(t,6) = value(V2(t));
                X(t,7) = value(p2(t));
                X(t,8) = value(u2(t));
                X(t,9) = value(s2(t));
                X(t,10)= q2;

                % POST-SOLVE CHECKS (ICC, t>1) 
                if (strcmp(bounds,"icc") || strcmp(bounds,"jcc-bon")) && t > 1
                    tol = 1e-6;
                    V1_t = X(t,1); V2_t = X(t,6);
                    assert(isfinite(V1_t) && isfinite(V2_t), 'V NaN/Inf at t=%d', t);
                    assert(V1_t >= s(1).min_V - tol && V1_t <= s(1).max_V + tol, ...
                           'V1 out of bounds at t=%d (%.4g)', t, V1_t);
                    assert(V2_t >= s(2).min_V - tol && V2_t <= s(2).max_V + tol, ...
                           'V2 out of bounds at t=%d (%.4g)', t, V2_t);
                end
    
            case "jcc-ssh"
                % Calculate Covariance Matrix
                rho = params.rho;
                Sigma_q = [std_hat(t,1)^2, rho*prod(std_hat(t,:));
                    rho*prod(std_hat(t,:)), std_hat(t,2)^2];

                if t == 1
                     % Standard solve
                    model = optimize(cons, -Objective, options);
        
                    X(t,1) = value(V1(t));
                    X(t,2) = value(p1(t));
                    X(t,3) = value(u1(t));
                    X(t,4) = value(s1(t));
                    X(t,5) = q1;
                    X(t,6) = value(V2(t));
                    X(t,7) = value(p2(t));
                    X(t,8) = value(u2(t));
                    X(t,9) = value(s2(t));
                    X(t,10)= q2;
                    
                else 
                    % Apply SSH Alg
                    vars = struct('V1',V1(t),'p1',p1(t),'u1',u1(t),'s1',s1(t), ...
                    'V2',V2(t),'p2',p2(t),'u2',u2(t),'s2',s2(t));

                    % Find Slater Point
                    q_mean = [q1; q2];
                    x_slater = findSlater(X(t-1,:), q_mean, s, c);
            
                    % Apply SSH Alg to determine topimal solution 
                    [cons, x_sol, phi_vals(t), alpha_vals(t,:)] = applySSH(cons, vars, t, ...
                                             X(t-1,:), q_mean, Sigma_q, ...
                                             x_slater, eps_t, s, Objective, options);
        
                    % Store result
                    X(t,1:10) = [x_sol(1) x_sol(2) x_sol(3) x_sol(4) q1 ...
                                 x_sol(5) x_sol(6) x_sol(7) x_sol(8) q2];
                end 
        end

    end

    obj = sum(X(:,2) + X(:,7));
    fprintf('Optimization complete.\n');

    % Simulation Report
    variable_report(framework, params.season, obj, X, N)

end


function variable_report(framework, season, obj, X, N)
    sd = 4;
    fprintf('\n');
    fprintf('--- VARIABLE REPORT ---\n');
    fprintf('Segments in PWL = %d\n', N);
    fprintf('Uncertainty Framework: %s\n', framework);
    fprintf('Season: %s\n', season);
    fprintf('\n');
    fprintf('Total Generation [MWh]: %g\n', round(obj, sd+1, 'significant'));
    fprintf('Generation 01 [MWh]: %g\n', round(sum(X(:,2)), sd+1, 'significant'));
    fprintf('Generation 02 [MWh]: %g\n', round(sum(X(:,7)), sd+1, 'significant'));
    fprintf('Generation Release 01 [m3]: %g\n', round(sum(X(:,3)), sd, 'significant'));
    fprintf('Generation Release 02 [m3]: %g\n', round(sum(X(:,8)), sd, 'significant'));
    fprintf('Spill Release 01 [m3]: %g\n', floor(round(sum(X(:,4)), sd-1, 'significant')));
    fprintf('Spill Release 02 [m3]: %g\n', floor(round(sum(X(:,9)), sd-1, 'significant')));
    fprintf('\n');
end


function h_ref = map_V_to_h(h, left_bounds, right_bounds, reference_values)
    N = length(left_bounds);
    if h <= left_bounds(1)
        h_ref = reference_values(1);
    elseif h >= right_bounds(end)
        h_ref = reference_values(N);
    else
        idx = find(left_bounds <= h, 1, 'last');
        idx = min(idx, N);
        h_ref = reference_values(idx);
    end
end

function q_norm = normalize_flow(flow, mean_val, std_val)
    q_norm = (flow - mean_val)/std_val;
end


function q = rescale_flow(flow_norm, mean_val, std_val)
    q = flow_norm*std_val + mean_val;
end


function [q1, q2, std_hat] = forecast_inflow(X, t, q1_s, lag, framework, params, s)
    
    std_hat = [0, 0];
    
    % Calculate previous release
    if t <= lag 
        X1_prev = s(1).min_ut; % min flow
    else 
        X1_prev = X(t-lag,3) + X(t-lag,4);
    end 
    
    %% Forecast Upstream Inflow
    switch framework
        case "det"
            % Calculate Q1
            q1 = q1_s(t + lag); % perfect foresight
            
            % Calculate Q2
            q2 = X1_prev; % lagged upstream release

        case "diu"
            % Forecast Q1
            [q1, std_hat(1)] = forecast_inflow_diu(q1_s(t), params); %q1_s is lagged historical inflow
            
            % Forecast Q2
            % (TEMP) will use q2_s(t) "local flow" in Columbia River
            [q2, std_hat(2)] = forecast_inflow_diu(X1_prev, params);
           
        
        case "ddu"
            if t <= lag
                % Forecast Q1
                [q1, std_hat(1)] = forecast_inflow_diu(q1_s(t), params); %q1_s is lagged historical inflow
                
                % Forecast Q2
                % (TEMP) will use q2_s(t) "local flow" in Columbia River
                [q2, std_hat(2)] = forecast_inflow_diu(X1_prev, params);
            else
                % Forecast Q1 (DIU - no upstream to condition on)
                [q1, std_hat(1)] = forecast_inflow_diu(q1_s(t), params);
               
                % Forecast Q2
                % (TEMP) will use q2_s(t) "local flow" for arg1
                [q2, std_hat(2)] = forecast_inflow_ddu(X1_prev, X(t-lag,10), X1_prev, params);
            end
    end
end

function [q_hat, std_hat_m3] = forecast_inflow_diu(q_prev, params)

    % (TEMP) scale to return forecasted std dev to original (dry)
    scale = 1; %3.14;
    
    % Normalize previous observed inflow
    q_prev_norm = normalize_flow(q_prev, params.inflow_mean, params.inflow_std);

    % Construct inflow forecast q_t = alpha_0 + alpha_1*q_{t-1}
    q_hat_norm = params.AR_const + params.AR_coef*q_prev_norm;

    % Rescale inflow forecast
    q_hat = rescale_flow(q_hat_norm, params.inflow_mean, params.inflow_std);

    % Rescale standard deviation back to inflow units
    std_hat_m3 = scale*params.AR_resid_std*params.inflow_std;
end


function [q_hat, std_hat_m3] = forecast_inflow_ddu(q_prev, q_pred_prev, outflow_prev, params)
    
    % Normalize predictors
    q_prev_norm = normalize_flow(q_prev, params.inflow_mean, params.inflow_std);
    outflow_prev_norm = normalize_flow(outflow_prev, params.outflow_mean, params.outflow_std);
    q_pred_prev_norm = normalize_flow(q_pred_prev, params.inflow_mean, params.inflow_std);

    % Calculate previous error term (norm)
    error_prev = abs(q_prev_norm - q_pred_prev_norm);
    norm_error_prev = (error_prev - params.error_mean)/params.error_std;

    % Forecast conditional variance using GARCH-X
    var_hat_norm = params.omega + params.alpha*(norm_error_prev^2) + params.gamma*outflow_prev_norm;
    std_hat_norm = sqrt(var_hat_norm);
    
    % Rescale standard deviation to norm inflow scale 
    std_hat = std_hat_norm*params.error_std + params.error_mean;

    % Construct inflow forecast
    q_hat_norm = params.constant + params.coef1*q_prev_norm + params.coef2*outflow_prev_norm;
    q_hat = rescale_flow(q_hat_norm, params.inflow_mean, params.inflow_std);

    % Rescale standard deviation back to inflow units 
    std_hat_m3 = std_hat*params.inflow_std;
end


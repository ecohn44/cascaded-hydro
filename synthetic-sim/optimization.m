function [sol, X, q, std_hat] = optimization(T, N, n, c, q_hist, lag, framework, modelparams, sysparams)

    % Initialize decision variable storage
    % X columns: 1=V1, 2=p1, 3=u1 
    %            4=V2, 5=p2, 6=u2, ...
    X = zeros(T,3*n);
    
    % Initialize streamflow vector 
    q = zeros(T,n);
    
    % Initialize process uncertainty vector 
    std_hat = zeros(T,1); 

    % Define decision variables (YALMIP)
    yalmip('clear');

    % Decision Variables Units 1,..., n
    V = sdpvar(T,n);
    p = sdpvar(T,n);
    u = sdpvar(T,n);

    % Objective
    cons = [];
    Objective = sum(sum(p));


    %% Build constraints
    for t = 1:T
        %% Forecast Upstream Inflow 
        switch framework
            case "det"
                q(t,1) = q_hist(t + lag); % perfect foresight
                std_hat(t) = 0;     % no variance
            
            case "diu"
                if t <= lag
                    q(t,1) = q_hist(t); % use prev inflow as predictor 
                    std_hat(t) = 0;
                else
                    [q(t,1), std_hat(t)] = forecast_inflow_diu(q_hist(t), modelparams); % DIU OlS model
                end
        end

        %% Estimate Downstream Inflow 
        %{
        for i = 2:n
            if t <= lag
                q(t,i) = 0; % no water has traveled yet 
            else
                 % inflow is lagged upstream release
                q(t,i) = X(t-lag,3*(i-2)+3); %+ rescale_flow(randn()*modelparams.AR_resid_std, modelparams.inflow_mean, modelparams.inflow_std); 
            end
        end 
        %}

        % Loop over all units
        for i = 1:n
            %% Static Constraints
            sp = sysparams(i);
        
            % Water Release Bounds
            cons = [cons, sp.min_ut <= u(t,i) <= sp.max_ut];
        
            % Feeder Capacity
            cons = [cons, 0 <= p(t,i) <= sp.F];
        
            % Volume Bounds
            cons = [cons, sp.min_V <= V(t,i) <= sp.max_V];

            %% Time-Varying Constraints
            if t == 1 % Initial conditions
                % Mass Balance
                %cons = [cons, V(t,i) == sp.V0]; 
                cons = [cons, V(t,i) + u(t,i) == sp.V0 + q(t,i)];
    
                % Ramp Rates (initialize at min release)
                cons = [cons, u(t,i) == sp.min_ut];
    
                % Power Production at t=1
                cons = [cons, p(t,i) == c * u(t,i) * sp.a * (sp.V0^sp.b)];
            
            else
                % Mass Balance
                cons = [cons, V(t,i) + u(t,i) == V(t-1,i) + q(t,i)];
    
                % Ramp Rate Up
                cons = [cons, u(t,i) <= sp.RR_up + u(t-1,i)];
    
                % Ramp Rate Down
                cons = [cons, u(t,i) >= u(t-1,i) - sp.RR_dn];
    
                % Power Production
                h_lookup = sp.a * X(t-1,3*(i-1)+1)^sp.b; % prev hydraulic head
                h_ref = map_V_to_h(h_lookup, sp.h_lbounds, sp.h_rbounds, sp.h_refvals);
                cons = [cons, p(t,i) == c * u(t,i) * h_ref];                
            end
        
        end

        % Solve optimization
        options = sdpsettings('solver','gurobi','verbose',0);
        sol = optimize(cons, -Objective, options);
        obj = 0;

         % Check solver status
        if sol.problem ~= 0
            fprintf('Solver returned problem = %d\n', sol.problem);
            disp(yalmiperror(sol.problem));
            % Optionally stop or break so you can inspect variables
            error('Solver failed at time t=%d (abort).', t);
        end
    
        % Store results for time t
        for i = 1:n
            X(t,3*(i-1)+1) = value(V(t,i));  % Volume
            X(t,3*(i-1)+2) = value(p(t,i));  % Power
            X(t,3*(i-1)+3) = value(u(t,i));  % Release
            obj = obj + value(p(t,i)); 
        end
    
    end
    fprintf('Optimization complete.\n');

    % Simulation Report
    variable_report(X, obj, n, N, framework, modelparams.season)

end


function variable_report(X, sol, n, N, framework, season)
    sd = 4;
    
    fprintf('\n');
    fprintf('--- VARIABLE REPORT ---\n');
    fprintf('Segments in PWL = %d\n', N);
    fprintf('Uncertainty Framework: %s\n', framework);
    fprintf('Season: %s\n', season);
    fprintf('\n');

    fprintf('Total Generation [MWh]: %g\n', round(sol, sd+1, 'significant'));

    for i = 1:n
        col_offset = 3*(i-1);  % 3 columns per unit: V, p, u

        gen_col = col_offset + 2;      % p
        release_col = col_offset + 3;  % u

        fprintf('Generation Unit %d [MWh]: %g\n', i, round(sum(X(:,gen_col)), sd+1, 'significant'));
        fprintf('Generation Release Unit %d [m3]: %g\n', i, round(sum(X(:,release_col)), sd, 'significant'));
    end

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

function [q_hat, e] = forecast_inflow_diu(q_prev, params)
    
    % Normalize previous observed inflow
    q_prev_norm = normalize_flow(q_prev, params.inflow_mean, params.inflow_std);
    
    % Sample and scale process uncertainty with historical std
    e = randn()*params.AR_resid_std;

    % Construct inflow forecast q_t = alpha_0 + alpha_1*q_{t-1} + e
    q_hat_norm = params.AR_const + params.AR_coef*q_prev_norm + e;

    % Rescale inflow forecast
    q_hat = rescale_flow(q_hat_norm, params.inflow_mean, params.inflow_std);
end


function [q_hat, std_hat] = forecast_inflow_ddu(q_prev, q_pred_prev, outflow_prev, params)
    
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
    % Rescale standard deviation 
    std_hat = std_hat_norm*params.error_std + params.error_mean;

    % Sample process uncertainty
    e = randn()*std_hat;

    % Construct inflow forecast
    q_hat_norm = params.constant + params.coef1*q_prev_norm + params.coef2*outflow_prev_norm + e;
    q_hat = rescale_flow(q_hat_norm, params.inflow_mean, params.inflow_std);
end


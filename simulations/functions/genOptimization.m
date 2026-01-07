%% Generalized n-unit Hydropower Optimization Model 
% =======================================================================
% Author: Eliza Cohn
% Description: Solves a single forward simulation/optimization loop over T time steps
% (non-anticipatory) with DIU/DDU/DET inflow models and JCC volume bounds.
% INPUTS
%   T           : Number of time periods in optimization horizon.
%   N           : Number of segments in the PWL head approximation (for reporting).
%   c           : Power conversion coefficient (unitless or W·s/m³ per your scaling).
%   q           : Exogenous inflow series (T + lag × n).
%   R           : Covariance (DIU) or correlation (DDU) matrix 
%   lag         : Lag length for inflow forecasting (typically 1).
%   scale       : Scaling factor for chance-constraint z-value.
%   framework   : Uncertainty type: "det", "diu", or "ddu".
%   bounds      : Volume bound type: "det" or "jcc-bon".
%   params      : Struct of uncertainty model parameters 
%                 (AR const/coef/std for DIU; omega/alpha/gamma for DDU).
%   s           : Struct array of system parameters for each unit i=1..n:
% 
% OUTPUTS
%   model       : Optimization model result for each time step.
%
%   obj         : Total generation over the horizon (sum of p_i).
%
%   X           : (T × 5n) matrix of optimal trajectories:
%                     For each unit i:
%                       col 5*(i-1)+1 = V_i(t)    Storage
%                       col 5*(i-1)+2 = p_i(t)    Power
%                       col 5*(i-1)+3 = u_i(t)    Turbine release
%                       col 5*(i-1)+4 = s_i(t)    Spill
%                       col 5*(i-1)+5 = q_i(t)    Predicted inflow
%
%   std_hat     : (T × n) matrix of predicted inflow standard deviations
%                 per unit under DIU/DDU frameworks.
%
%   V_eff       : (T × 2n) matrix of effective volume bounds with chance shifts:
%                     For each unit i:
%                        col 2*i-1 = max_V_eff(i,t)
%                        col 2*i   = min_V_eff(i,t)
% -------------------------------------------------------------------------


function [model, obj, X, std_hat, V_eff, phi_vals, alpha_vals] = genOptimization(T, N, c, q, lag, scale, framework, bounds, params, s)


    % Number of units
    n = numel(s);
    phi_vals   = zeros(T,1);
    alpha_vals = zeros(T,n);

    % Indexing scheme:
    % For unit i, columns in X are:
    %   idxV = 5*(i-1)+1   -> V_i
    %   idxP = 5*(i-1)+2   -> p_i
    %   idxU = 5*(i-1)+3   -> u_i
    %   idxS = 5*(i-1)+4   -> s_i (spill)
    %   idxQ = 5*(i-1)+5   -> q_i (predicted inflow)

    % Initialize storage
    X          = zeros(T, 5*n);
    std_hat    = zeros(T, n);
    std_safety = zeros(T, n);
    V_eff      = zeros(T, 2*n);  % [max_V_eff(i), min_V_eff(i)]

    % Risk level
    eps_t = 0.05;

    % YALMIP reset
    yalmip('clear');

    
    lambda = zeros(n,1);
    for i = 1:n
        % Use V0 as the reference point for the marginal value of head
        v_target = 0.05 * s(i).max_V;
        
        % Gradient: dh/dV
        dh_dV = s(i).a * s(i).b * (v_target ^ (s(i).b - 1));
        
        % Price: c * Expected_Release * Gradient
        u_expected = s(i).max_ut; 

        tie_breaker = 1.0 + (0.1 * (n - i + 1));
        
        % Scale down slightly to ensure Power >> Head Value (prioritize generation)
        lambda(i) = (c * u_expected * dh_dV)*tie_breaker ; 
    end
    


    %% Non-Anticipatory Optimization Framework 
    for t = 1:T
        
        %% Forecast Inflow for all units
        % q_t(i) = predicted mean inflow into unit i
        % std_hat(t,i) = std dev estimate
        % std_safety(t,i) = safety-scaled std (for bounds)
        [q_t, std_hat(t,:), std_safety(t,:)] = ...
            forecast_inflow(X, t, q, lag, framework, params, s);

        V_min_shift = zeros(1, n);
        V_max_shift = zeros(1, n);

        %% Determine Volume Bounds (chance constraints)
        switch bounds
            case {"det"}
                % Deterministic / SSH: no shifts

            case {"jcc-bon"}
                % Bonferroni z-score
                z = scale * norminv(1 - (eps_t/n));

                % Use safety-scaled std for all units
                V_min_shift =  z .* std_safety(t,:);
                V_max_shift = -z .* std_safety(t,:);

                % Store effective bounds and apply shift in constraints
                for i = 1:n
                    V_eff(t,2*i-1) = s(i).max_V + V_max_shift(i);
                    V_eff(t,2*i)   = s(i).min_V + V_min_shift(i);
                end

            case {"jcc-ssh"}
                % SSH z-score
                z = scale * norminv(1 - (eps_t));

                % Store effective vounds but don't apply shift directly
                for i = 1:n
                    V_eff(t,2*i-1) = s(i).max_V - z*std_safety(t,i);
                    V_eff(t,2*i)   = s(i).min_V + z*std_safety(t,i);
                end
        end

        %% Define Decision Variables (for this time step)
        V  = sdpvar(n,1);
        p  = sdpvar(n,1);
        u  = sdpvar(n,1);
        sp = sdpvar(n,1);  % spill

        cons = [];

        %{
        lambda = zeros(n,1);

        % Calculate marginal cost of volume 
        for i = 1:n
            
            % Retrieve previous volume 
            if t == 1
                v_current = s(i).V0;
            else
                idxV = 5*(i-1)+1;
                v_current = max(X(t-1, idxV), 0.01); 
            end
            
            % Calculate gradient dh/dV = a * b * V^(b-1)
            dh_dV = s(i).a * s(i).b * (v_current ^ (s(i).b - 1));
            
            % Calculate price of head
            u_expected = s(i).max_ut * 0.5; 
            lambda(i)  = c * u_expected * dh_dV;
        end
        %}

        % Maximizing (Power - Spill) + (Value of Stored Head)
        Objective = sum(p) - sum(sp) + sum(lambda .* V);  

        %% Static Constraints
        for i = 1:n
            % Water Release Bounds
            cons = [cons, s(i).min_ut <= u(i) <= s(i).max_ut];

            % Feeder Capacity
            cons = [cons, 0 <= p(i) <= s(i).F];

            % Spill Bounds
            cons = [cons, 0 <= sp(i)];
        end

        %% Time-Varying Constraints
        if t == 1
            
            % Initial conditions
            for i = 1:n
                % Mass balance: V_i + u_i + s_i = V0_i + q_t(i)
                cons = [cons, V(i) + u(i) + sp(i) == s(i).V0 + q_t(i)];

                % Initial ramp: fix at min_ut
                cons = [cons, u(i) == s(i).min_ut];

                % Power production (head at initial volume)
                h0 = s(i).a * (s(i).V0^s(i).b);
                cons = [cons, p(i) == c*u(i)*h0];

                % Volume bounds (with chance-constrained shift)
                cons = [cons, s(i).min_V + V_min_shift(i) <= s(i).V0 + q_t(i) - (u(i) + sp(i)) ...
                    <= s(i).max_V + V_max_shift(i)];
            end

        else
            % t > 1: depend on previous X
            for i = 1:n
                idxV = 5*(i-1)+1;
                idxU = 5*(i-1)+3;

                V_prev = X(t-1, idxV);
                u_prev = X(t-1, idxU);

                % Mass balance
                cons = [cons, V(i) + u(i) + sp(i) == V_prev + q_t(i)];

                % Ramp rates
                cons = [cons, s(i).RR_dn <= u(i) - u_prev <= s(i).RR_up];

                % Power production using previous volume's head (with PWL mapping)
                h_lookup = s(i).a * V_prev^s(i).b;
                h_ref    = map_V_to_h(h_lookup, s(i).h_lbounds, s(i).h_rbounds, s(i).h_refvals);
                cons = [cons, p(i) == c*u(i)*h_ref];

                % Enforce volume bounds with chance-constrained shift 
                cons = [cons, s(i).min_V + V_min_shift(i) <= V_prev + q_t(i) - (u(i) + sp(i)) ...
                                    <= s(i).max_V + V_max_shift(i)];
                
            end
        end

        %% Solve optimization at time t
        options = sdpsettings('solver','gurobi','verbose',0, 'gurobi.Seed', 1, 'gurobi.Threads', 1);
        
        % Use Gurobi for standard solve 
        switch bounds
            case {"det","icc","jcc-bon"}
                model = optimize(cons, -Objective, options);
        
                % Abort early on infeasible/numerical issues
                if model.problem ~= 0
                    fprintf('\n*** Solver issue at t = %d ***\n', t);
                    fprintf('  problem code : %d\n', model.problem);
                    fprintf('  yalmiperror  : %s\n', yalmiperror(model.problem));
                    fprintf('  model.info   : %s\n', model.info);
                    error('Stopping after first solver error at t=%d.', t);
                end
        
                % Store results for time t
                for i = 1:n
                    idxV = 5*(i-1)+1;
                    idxP = 5*(i-1)+2;
                    idxU = 5*(i-1)+3;
                    idxS = 5*(i-1)+4;
                    idxQ = 5*(i-1)+5;
        
                    X(t,idxV) = value(V(i));
                    X(t,idxP) = value(p(i));
                    X(t,idxU) = value(u(i));
                    X(t,idxS) = value(sp(i));
                    X(t,idxQ) = q_t(i);   % predicted inflow
                end

            case {"jcc-ssh"}

                % Calculate Covariance Matrix
                switch framework 
                    case {"det"}
                        % No inflow uncertainty: zero covariance
                        Sigma_q = zeros(n,n);
                
                    case {"diu"}
                        % DIU: use static covariance estimated offline
                        Sigma_q = params.Sigma_diu;          % n×n covariance
                
                    case {"ddu"}
                        % DDU: time-varying covariance Sigma_t = D_t Rcorr D_t
                        sigma_t = std_safety(t,:).';         % n×1 vector of sigma_i,t
                        D_t     = diag(sigma_t);             % n×n
                        Sigma_q = D_t * params.Rcorr * D_t;  % n×n time-varying Sigma_t
                end

                if t == 1
                    % Standard solve
                    model = optimize(cons, -Objective, options);
        
                    % Store results for time t
                    fprintf('@t=%d \n', t);
                    for i = 1:n
                        idxV = 5*(i-1)+1;
                        idxP = 5*(i-1)+2;
                        idxU = 5*(i-1)+3;
                        idxS = 5*(i-1)+4;
                        idxQ = 5*(i-1)+5;
            
                        X(t,idxV) = value(V(i));
                        X(t,idxP) = value(p(i));
                        X(t,idxU) = value(u(i));
                        X(t,idxS) = value(sp(i));
                        X(t,idxQ) = q_t(i);   % predicted inflow

                        fprintf('    Unit %d: \n', i);
                        fprintf('       V_t=%.4f, u_t=%.4f, p_t=%.4f, s_t=%.4f,  q_t=%.4f\n', ...
                            X(t,idxV), X(t,idxU), X(t,idxP), X(t,idxS), q_t(i));
                    end
                    
                else 
                    % Apply SSH Alg
                    vars = struct( ...
                        'V',  V,  ...   % n×1 sdpvar of volumes at time t
                        'p',  p,  ...   % n×1 power
                        'u',  u,  ...   % n×1 turbine releases
                        's',  sp ...    % n×1 spills
                    );

                    % Find Slater Point
                    fprintf('\n@t=%d SSH Initialization: \n', t)
                    x_slater = findSlater(X(t-1,:), q_t, s, c, V_eff(t,:));
            
                    target_phi = (1 - eps_t);
                
                    % Apply SSH  
                    [cons, x_sol, phi_k, alpha_k] = applySSH(cons, vars, t, ...
                        X(t-1,:), q_t, Sigma_q, x_slater, target_phi, s, ...
                        Objective, options);
                
                    % Store latest values for this time step
                    phi_vals(t)      = phi_k;
                    alpha_vals(t, :) = alpha_k;
                    viol       = max(0, target_phi - phi_k);    
        
                    % Store result into X(t,:)
                    for i = 1:n
                        idxV = 5*(i-1)+1;
                        idxP = 5*(i-1)+2;
                        idxU = 5*(i-1)+3;
                        idxS = 5*(i-1)+4;
                        idxQ = 5*(i-1)+5;
                    
                        base = 4*(i-1);       % position in x_sol
                        X(t,idxV) = x_sol(base+1);
                        X(t,idxP) = x_sol(base+2);
                        X(t,idxU) = x_sol(base+3);
                        X(t,idxS) = x_sol(base+4);
                        X(t,idxQ) = q_t(i);   % same q_t you just used

                        fprintf('    Unit %d: ', i);
                        fprintf('       V_t=%.4f, u_t=%.4f, p_t=%.4f, s_t=%.4f,  q_t=%.4f\n', ...
                            X(t,idxV), X(t,idxU), X(t,idxP), X(t,idxS), q_t(i));
                    end
                end 
        end 
    end 

    %% Objective: Total power production over all units
    obj = 0;
    for i = 1:n
        idxP = 5*(i-1)+2;
        obj = obj + sum(X(:,idxP));
    end

    fprintf('Optimization complete.\n');

    % Simulation Report
    variable_report(framework, params.season, obj, X, N, s)

end


function variable_report(framework, season, obj, X, N, s)
    sd = 4;
    n  = numel(s);

    total_gen = 0;
    for i = 1:n
        idxP = 5*(i-1)+2;
        total_gen = total_gen + sum(X(:,idxP));
    end

    fprintf('\n');
    fprintf('--- VARIABLE REPORT ---\n');
    fprintf('Segments in PWL = %d\n', N);
    fprintf('Uncertainty Framework: %s\n', framework);
    fprintf('Season: %s\n', season);
    fprintf('\n');
    fprintf('Total Generation [MWh]: %g\n', round(total_gen, sd+1, 'significant'));

    for i = 1:n
        idxP = 5*(i-1)+2;
        idxU = 5*(i-1)+3;
        idxS = 5*(i-1)+4;

        fprintf('Generation %02d [MWh]: %g\n', i, ...
            round(sum(X(:,idxP)), sd+1, 'significant'));
        fprintf('Generation Release %02d [m3]: %g\n', i, ...
            round(sum(X(:,idxU)), sd, 'significant'));
        fprintf('Spill Release %02d [m3]: %g\n', i, ...
            floor(round(sum(X(:,idxS)), sd-1, 'significant')));
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


function [q_t, std_hat, std_safety] = forecast_inflow(X, t, q, lag, framework, params, s)
    
    n = numel(s);

    q_t        = zeros(1,n);
    std_hat    = zeros(1,n);
    std_safety = zeros(1,n);

    % Compute previous outflows (t-lag) for each unit
    outflow_prev = zeros(1,n);
    
    if t <= lag
        % Fallback: minimum turbine release, zero spill
        for i = 1:n
            outflow_prev(i) = s(i).min_ut;
        end
    else
        row = t - lag;
        for i = 1:n
            idxU = 5*(i-1) + 3;   % u index in X
            idxS = 5*(i-1) + 4;   % s index in X
            outflow_prev(i) = X(row, idxU) + X(row, idxS);
        end
    end

    % Calculate inflow for each unit
    for i = 1:n
        
        % Local inflow at time t for unit i (already lagged)
        q_loc_prev = q(t, i);         
        
        % Upstream lagged outflow for cascade in series
        if i > 1
            up_out_prev = outflow_prev(i-1);
        else
            up_out_prev = outflow_prev(1);   
        end

        switch framework

            case "det"
                % Deterministic: just use the DIU mean, ignore std
                [q_t(i), ~] = forecast_inflow_diu(i, q_loc_prev, params);
                % std_hat and std_safety remain zero

            case "diu"
                % DIU: each unit uses DIU forecast and same std for safety
                [q_t(i), std_hat(i)] = forecast_inflow_diu(i, q_loc_prev, params);
                std_safety(i) = std_hat(i);

            case "ddu"
                if i == 1
                    
                    % Unit 1: no upstream operator to condition on (DIU)
                    [q_t(i), std_hat(i)] = forecast_inflow_diu(i, q_loc_prev, params);
                    std_safety(i) = std_hat(i);

                else

                    if t <= lag
                        % GARCH initialization
                        q_pred_prev_norm = q_loc_prev;
                    else
                        % Pull previous predicted inflow from X 
                        idxQ         = 5*(i-1) + 5;
                        q_pred_prev  = X(t-lag, idxQ);   % t-1, not t-lag, for GARCH error
                        q_pred_prev_norm = q_pred_prev;
                    end

                    % DDU variance model 
                    q_prev_norm       = q_loc_prev;      % previous observed inflow
                    outflow_prev_norm = up_out_prev;     % lagged upstream outflow

                    [q_t(i), std_hat(i), std_safety(i)] = ...
                        forecast_inflow_ddu(q_prev_norm, q_pred_prev_norm, ...
                                            outflow_prev_norm, params);
                end

            otherwise
                error('Unknown framework: %s', framework);
        end
    end
end


function [q_hat, std_hat] = forecast_inflow_diu(i, q_prev, params)
    
    % Construct inflow forecast q_t = alpha_0 + alpha_1*q_{t-1}
    q_hat = params.AR_const + params.AR_coef*q_prev;

    % Unit standard deviation
    std_hat = params.sigma_diu(i);
end


% (TEMP) add back in normalization for full river 
function [q_hat, std_hat_estimate, std_hat_safety] = forecast_inflow_ddu(q_prev_norm, q_pred_prev_norm, outflow_prev_norm, params)

    % Sensitivity to upstream release
    scale_exog = 1;

    % Calculate previous error term
    norm_error_prev = abs(q_prev_norm - q_pred_prev_norm);
    
    % Forecast conditional variance using GARCH-X
    var_hat_norm     = params.omega + params.alpha*(norm_error_prev^2) + params.gamma*outflow_prev_norm;
    std_hat_estimate = sqrt(var_hat_norm);

    % Safety-bound estimate of conditional variance using GARCH-X
    var_hat_norm     = params.omega + params.alpha*(norm_error_prev^2) + scale_exog*params.gamma*outflow_prev_norm;
    std_hat_safety   = sqrt(var_hat_norm);
    
    % Construct inflow forecast from OLS model
    q_hat = params.constant + params.coef1*q_prev_norm + params.coef2*outflow_prev_norm;
  
end
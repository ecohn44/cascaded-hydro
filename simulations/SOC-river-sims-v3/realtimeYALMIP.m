function [result, obj, X, std_hat, phi_val, alpha_vals] = realtimeYALMIP(t, c, kappa, eps, I_prev, q_error, V_prev, u_prev, V_ref, theta, lag, up_release, sys, model, bounds, framework, tracking)
% =========================================================================
% HEAD FORMULA (linear):
%   h_i = a_i * V_i + b_i          (sys(i).a = slope, sys(i).b = intercept)
%
% POWER FORMULA:
%   p_i = c * h_i * u_i
%       = c * (a_i * V_i + b_i) * u_i
%       = c * a_i * w_i + c * b_i * u_i    where w_i = V_i * u_i
%
% =========================================================================

    % Clear YALMIP internal variable database
    yalmip('clear');

    n = numel(sys);

    Rcorr = [ 1.000000, -0.033090,  0.020604;
             -0.033090,  1.000000,  0.014311;
              0.020604,  0.014311,  1.000000]; 

    % Step 1: Forecast inflow and estimate standard deviation
    q_t = I_prev(:);
    if t > lag
        for i = 2:n
            j = i - 1;  % model 1=JDA, 2=TDA, 3=BON
            q_t(i) = model(j).coef0 + model(j).coef1 * I_prev(i) + model(j).coef2 * up_release(i-1);
        end
    end
    std_hat = forecast_error(t, kappa, q_error, up_release, framework, model, sys);

    % Step 2: Compute effective volume bounds (JCC volume shift)
    switch bounds
        case {"det", "jcc-ssh"}
            for i = 1:n
                sys(i).V_eff_max = sys(i).max_V;
                sys(i).V_eff_min = sys(i).min_V;
            end
        case {"jcc-bon"}

            z_score     =  norminv(1 - (eps / (2 * n)));
            std_V   = [sys.kV]' .* std_hat(:); % convert to volume

            for i = 1:n
                sys(i).V_eff_min = sys(i).min_V + z_score*std_V(i);
                sys(i).V_eff_max = sys(i).max_V;
            end
    end

    % Step 3: Declare YALMIP decision variables
    V  = sdpvar(n, 1);
    p  = sdpvar(n, 1);
    u  = sdpvar(n, 1);
    sp = sdpvar(n, 1);
    d  = sdpvar(n, 1);
    w  = sdpvar(n, 1);

    % Step 4: Build constraint set
    cons = buildConstraints(t, n, V, p, u, sp, d, w, V_prev, u_prev, q_t, V_ref, c, sys);

    % Step 5: Build objective function
    P_base = sum([sys.F]);
    Objective   = -sum(p)/P_base + 1e-4 * sum(sp) + sum(theta .* d);

    % Step 6: Solve
    phi_val    = NaN;
    alpha_vals = zeros(n, 1);

    ops = sdpsettings('solver','gurobi','verbose',0, ...
                      'gurobi.TimeLimit',10,'gurobi.Seed',1,'gurobi.Threads',1);

    switch bounds
        case {"det", "jcc-bon"}
            [result, obj, X] = solveAndExtract(cons, Objective, V, p, u, sp, q_t, ops);

        case "jcc-ssh"
            % Build Sigma_q
            switch framework
                case "det",  Sigma_q = zeros(n,n);
                case "diu",  Sigma_q = model.AR_coef*eye(n);
                case "ddu"
                    D_t     = diag(std_hat(:));
                    Sigma_q = D_t * Rcorr * D_t;
            end

            if t == 1
                [result, obj, X] = solveAndExtract(cons, Objective, V, p, u, sp, q_t, ops);
            else
                % Construct X_prev row (1 x 5n) expected by applySSH
                X_prev_row = zeros(1, 5*n);
                for i = 1:n
                    X_prev_row(5*(i-1)+1) = V_prev(i);
                    X_prev_row(5*(i-1)+3) = u_prev(i);
                end

                vars     = struct('V', V, 'p', p, 'u', u, 's', sp);
                x_slater = findSlater(X_prev_row, q_t, sys, c, model.season);

                [~, x_sol, phi_val, alpha_vals, ~] = applySSH(cons, vars, t, ...
                    X_prev_row, q_t', Sigma_q, x_slater, 1 - eps, sys, Objective, ops);

                % Pack outputs from x_sol
                result = struct('problem', 0);
                X      = zeros(n, 5);
                for i = 1:n
                    base   = 4*(i-1);
                    X(i,:) = [x_sol(base+1), x_sol(base+2), x_sol(base+3), x_sol(base+4), q_t(i)];
                end
                obj = sum(X(:,2));
            end
    end

end



function cons = buildConstraints(t, n, V, p, u, sp, d, w, V_prev, u_prev, q_t, V_ref, c, sys)

    cons = [];

    V_lower = V_ref(1:n);
    V_upper = V_ref(n+1:2*n);

    for i = 1:n

        VL = sys(i).V_eff_min;
        VU = sys(i).V_eff_max;
        uL = sys(i).min_ut;
        uU = sys(i).max_ut;

        % (B1): Storage bounds
        cons = [cons, VL <= V(i) <= VU];

        % (B2): Power bounds
        cons = [cons, 0 <= p(i) <= sys(i).F];

        % (B3) Turbine release bounds + ramp rates
        if t == 1
            cons = [cons, u(i) == uL];
        else
            cons = [cons, uL <= u(i) <= uU];
            cons = [cons, u(i) >= u_prev(i) + sys(i).RR_dn];
            cons = [cons, u(i) <= u_prev(i) + sys(i).RR_up];
        end

        % (B4-B6) Non-negativity and McCormick auxiliary bounds
        cons = [cons, sp(i) >= 0, d(i) >= 0, 0 <= w(i) <= VU*uU];

        % (C1) Mass balance
        cons = [cons, V(i) == V_prev(i) + sys(i).kV*(q_t(i) - u(i) - sp(i))];

        % (C4) McCormick envelope for w(i) = V(i) * u(i)
        cons = [cons, ...
            w(i) >= VL*u(i) + uL*V(i) - VL*uL, ...    % MC1
            w(i) >= VU*u(i) + uU*V(i) - VU*uU, ...    % MC2
            w(i) <= VU*u(i) + uL*V(i) - VU*uL, ...    % MC3
            w(i) <= VL*u(i) + uU*V(i) - VL*uU  ...    % MC4
        ];

        % (C5) Power production
        cons = [cons, p(i) == c * sys(i).a * w(i) + c * sys(i).b * u(i)];

        % (C6-C7) Tracking error
        cons = [cons, V(i) - d(i) <= V_upper(i), -V(i) - d(i) <= -V_lower(i)];

    end
end


function [result, obj, X] = solveAndExtract(Constraints, Objective, V, p, u, sp, q_t, ops)

    result   = optimize(Constraints, Objective, ops);
    feasible = ismember(result.problem, [0, 3, 5]);

    obj = NaN;
    X   = zeros(numel(q_t), 5);
    if feasible
        X   = [value(V), value(p), value(u), value(sp), q_t(:)];
        obj = sum(value(p));
    end
end

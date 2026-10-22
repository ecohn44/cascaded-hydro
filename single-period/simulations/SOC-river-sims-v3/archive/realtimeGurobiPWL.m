function [result, obj, X, std_hat] = realtimeGurobiPWL(t, c, kappa, eps, I_prev, q_error, V_prev, u_prev, V_ref, theta, lag, up_release, sys, model, bounds, framework, tracking)
% =========================================================================
% HEAD FORMULA (linear):
%   h_i = a_i * V_i + b_i        (sys(i).a = slope, sys(i).b = intercept)
%
% POWER FORMULA:
%   p_i = c * h_i * u_i = c*a_i*w_i + c*b_i*u_i   where w_i = V_i * u_i
%
% McCormick envelope of w_i = V_i * u_i over [VL, VU] x [uL, uU]:
%   MC1:  w >= VL*u + uL*V - VL*uL
%   MC2:  w >= VU*u + uU*V - VU*uU
%   MC3:  w <= VU*u + uL*V - VU*uL
%   MC4:  w <= VL*u + uU*V - VL*uU
% =========================================================================

    n = numel(sys);

    % 1: Forecast inflow and estimate error
    q_t = I_prev(:);
    if t > lag
        for i = 2:n
            q_t(i) = model.coef0 + model.coef1 * I_prev(i) + model.coef2 * up_release(i-1);
        end
    end
    std_hat = forecast_error(t, kappa, q_error, up_release, framework, model, sys);

    % 2: Calculate volume shift
    switch bounds
        case {"det", "jcc-ssh"}
            for i = 1:n
                sys(i).V_eff_max = sys(i).max_V;
                sys(i).V_eff_min = sys(i).min_V;
            end
        case {"jcc-bon"}
            z_score     =  norminv(1 - (eps / (2*n)));
            V_min_shift =  z_score .* std_hat;
            V_max_shift = -z_score .* std_hat;
            for i = 1:n
                sys(i).V_eff_max = sys(i).max_V + V_max_shift(i);
                sys(i).V_eff_min = max(0, sys(i).min_V + V_min_shift(i));
            end
    end

    % 3: Decision variable layout
    nVars  = 6 * n;
    idx_V  = @(i)  (i - 1);
    idx_p  = @(i)  n   + (i - 1);
    idx_u  = @(i)  2*n + (i - 1);
    idx_sp = @(i)  3*n + (i - 1);
    idx_d  = @(i)  4*n + (i - 1);
    idx_w  = @(i)  5*n + (i - 1);   % w(i) = V(i) * u(i)  [McCormick auxiliary]

    % 4: Variable bounds and objective coefficients
    [lb, ub, obj_coeff] = buildVariables(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_w, theta, sys);

    % 5: Linear constraints
    [A, rhs, sense] = buildLinearConstraints(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_w, V_prev, u_prev, q_t, V_ref, c, sys);

    % 6: Solve
    result = solveModel(nVars, obj_coeff, A, rhs, sense, lb, ub, t);

    % 7: Extract solution
    [obj, X] = extractSolution(result, t, n, idx_V, idx_p, idx_u, idx_sp, q_t, V_ref, sys, tracking);

end


% =========================================================================
% Helper: variable bounds and objective coefficients
% =========================================================================
function [lb, ub, obj_coeff] = buildVariables(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_w, theta, sys)

    lb        = zeros(nVars, 1);
    ub        =  inf(nVars, 1);
    obj_coeff =  zeros(nVars, 1);

    for i = 1:n

        % B1: Storage V(i)
        lb(idx_V(i)+1)  = sys(i).V_eff_min;
        ub(idx_V(i)+1)  = sys(i).V_eff_max;

        % B2: Power output p(i)
        lb(idx_p(i)+1)  = 0;
        ub(idx_p(i)+1)  = sys(i).F;
        obj_coeff(idx_p(i)+1) = -1;

        % B3: Turbine release u(i)
        if t == 1
            lb(idx_u(i)+1) = sys(i).min_ut;
            ub(idx_u(i)+1) = sys(i).min_ut;
        else
            lb(idx_u(i)+1) = sys(i).min_ut;
            ub(idx_u(i)+1) = sys(i).max_ut;
        end

        % B4: Spill sp(i)
        lb(idx_sp(i)+1) = 0;
        ub(idx_sp(i)+1) = inf;
        obj_coeff(idx_sp(i)+1) = 1e-4;

        % B5: Tracking error d(i)
        lb(idx_d(i)+1)  = 0;
        ub(idx_d(i)+1)  = inf;
        obj_coeff(idx_d(i)+1) = theta / (sys(i).max_V - sys(i).min_V);

        % B6: McCormick auxiliary w(i) = V(i)*u(i) in [VL*uL, VU*uU]
        lb(idx_w(i)+1)  = sys(i).V_eff_min * sys(i).min_ut;
        ub(idx_w(i)+1)  = sys(i).V_eff_max * sys(i).max_ut;

    end
end


% =========================================================================
% Helper: linear constraints
% =========================================================================
function [A, rhs, sense] = buildLinearConstraints(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_w, V_prev, u_prev, q_t, V_ref, c, sys)

    col = @(x) x(:);

    rows  = zeros(0,1);
    cols  = zeros(0,1);
    vals  = zeros(0,1);
    rhs   = zeros(0,1);
    sense = '';
    row   = 0;

    for i = 1:n

        VL = sys(i).V_eff_min;
        VU = sys(i).V_eff_max;
        uL = sys(i).min_ut;
        uU = sys(i).max_ut;

        % (C1) Mass balance: V(i) + u(i) + sp(i) = V_prev(i) + q_t(i)
        row = row + 1;
        rows  = [rows;  row;         row;         row        ];
        cols  = [cols;  idx_V(i)+1;  idx_u(i)+1;  idx_sp(i)+1];
        vals  = [vals;  1;           1;           1          ];
        rhs   = [rhs;   V_prev(i) + q_t(i)                   ];
        sense = [sense; '='                                   ];

        if t > 1
            % (C2) Ramp-rate lower bound: u(i) >= u_prev(i) + RR_dn
            row = row + 1;
            rows  = [rows;  row        ];
            cols  = [cols;  idx_u(i)+1 ];
            vals  = [vals;  -1         ];
            rhs   = [rhs;   -(u_prev(i) + sys(i).RR_dn)];
            sense = [sense; '<'        ];

            % (C3) Ramp-rate upper bound: u(i) <= u_prev(i) + RR_up
            row = row + 1;
            rows  = [rows;  row        ];
            cols  = [cols;  idx_u(i)+1 ];
            vals  = [vals;  1          ];
            rhs   = [rhs;   u_prev(i) + sys(i).RR_up];
            sense = [sense; '<'        ];
        end

        % (C4) McCormick envelope for w(i) = V(i) * u(i)
        % MC1: -w + VL*u + uL*V <= VL*uL
        row = row + 1;
        rows  = [rows;  row;         row;         row        ];
        cols  = [cols;  idx_w(i)+1;  idx_u(i)+1;  idx_V(i)+1 ];
        vals  = [vals;  -1;          VL;          uL         ];
        rhs   = [rhs;   VL * uL                               ];
        sense = [sense; '<'                                   ];

        % MC2: -w + VU*u + uU*V <= VU*uU
        row = row + 1;
        rows  = [rows;  row;         row;         row        ];
        cols  = [cols;  idx_w(i)+1;  idx_u(i)+1;  idx_V(i)+1 ];
        vals  = [vals;  -1;          VU;          uU         ];
        rhs   = [rhs;   VU * uU                               ];
        sense = [sense; '<'                                   ];

        % MC3: w - VU*u - uL*V <= -VU*uL
        row = row + 1;
        rows  = [rows;  row;         row;         row        ];
        cols  = [cols;  idx_w(i)+1;  idx_u(i)+1;  idx_V(i)+1 ];
        vals  = [vals;  1;          -VU;         -uL         ];
        rhs   = [rhs;  -VU * uL                               ];
        sense = [sense; '<'                                   ];

        % MC4: w - VL*u - uU*V <= -VL*uU
        row = row + 1;
        rows  = [rows;  row;         row;         row        ];
        cols  = [cols;  idx_w(i)+1;  idx_u(i)+1;  idx_V(i)+1 ];
        vals  = [vals;  1;          -VL;         -uU         ];
        rhs   = [rhs;  -VL * uU                               ];
        sense = [sense; '<'                                   ];

        % (C5) Power production: p(i) = c*a_i*w(i) + c*b_i*u(i)
        row = row + 1;
        rows  = [rows;  row;         row;                  row        ];
        cols  = [cols;  idx_p(i)+1;  idx_w(i)+1;           idx_u(i)+1 ];
        vals  = [vals;  1;          -c * sys(i).a;         -c * sys(i).b];
        rhs   = [rhs;   0                                               ];
        sense = [sense; '='                                             ];

        % (C6) Tracking error upper: V(i) - d(i) <= V_upper(i)
        V_upper = V_ref(n+1:2*n);
        row = row + 1;
        rows  = [rows;  row;         row        ];
        cols  = [cols;  idx_V(i)+1;  idx_d(i)+1 ];
        vals  = [vals;  1;          -1          ];
        rhs   = [rhs;   V_upper(i)               ];
        sense = [sense; '<'                      ];

        % (C7) Tracking error lower: -V(i) - d(i) <= -V_lower(i)
        V_lower = V_ref(1:n);
        row = row + 1;
        rows  = [rows;  row;         row        ];
        cols  = [cols;  idx_V(i)+1;  idx_d(i)+1 ];
        vals  = [vals;  -1;         -1          ];
        rhs   = [rhs;  -V_lower(i)               ];
        sense = [sense; '<'                      ];

    end

    nCons = row;
    A = sparse(col(rows), col(cols), col(vals), nCons, nVars);
end


% =========================================================================
% Helper: solve model (pure LP — no genconpow, no FuncNonlinear)
% =========================================================================
function result = solveModel(nVars, obj_coeff, A, rhs, sense, lb, ub, t_period)

    grb_model.modelname  = sprintf('hydroRT_t%d', t_period);
    grb_model.modelsense = 'min';
    grb_model.obj        = obj_coeff(:);
    grb_model.A          = A;
    grb_model.rhs        = rhs(:);
    grb_model.sense      = sense;
    grb_model.lb         = lb(:);
    grb_model.ub         = ub(:);
    grb_model.vtype      = repmat('C', nVars, 1);

    params.OutputFlag = 0;
    params.Seed       = 1;
    params.TimeLimit  = 10;
    params.Threads    = 1;

    result = gurobi(grb_model, params);

    if strcmp(result.status, 'INFEASIBLE')
        iis = gurobi_iis(grb_model);
        iis_rows = find(iis.Arows);
        for r = iis_rows'
            fprintf('Row %d | sense=%s | rhs=%.6f\n', r, grb_model.sense(r), grb_model.rhs(r));
            [~, icols, ivals] = find(grb_model.A(r,:));
            for k = 1:length(icols)
                fprintf('  col %d  coeff=%.6f  lb=%.6f  ub=%.6f\n', ...
                    icols(k), ivals(k), grb_model.lb(icols(k)), grb_model.ub(icols(k)));
            end
        end
    end

end


% =========================================================================
% Helper: extract solution
% =========================================================================
function [obj, X] = extractSolution(result, t, n, idx_V, idx_p, idx_u, idx_sp, q_t, V_ref, sys, tracking)

    obj = NaN;
    X   = zeros(n, 5);

    has_solution = ismember(result.status, {'OPTIMAL','SUBOPTIMAL','TIME_LIMIT'}) ...
                   && isfield(result,'x') && ~isempty(result.x);

    if has_solution
        x = result.x;

        V_out  = zeros(n,1);
        p_out  = zeros(n,1);
        u_out  = zeros(n,1);
        sp_out = zeros(n,1);

        for i = 1:n
            V_out(i)  = x(idx_V(i)+1);
            p_out(i)  = x(idx_p(i)+1);
            u_out(i)  = x(idx_u(i)+1);
            sp_out(i) = x(idx_sp(i)+1);
        end

        % Tracking error (diagnostics only)
        track_err = zeros(n,1);
        switch tracking
            case "mean"
                for i = 1:n
                    vol_range     = sys(i).max_V - sys(i).min_V;
                    track_err(i)  = abs(V_out(i) - V_ref(i)) / vol_range;
                end
            case "envelope"
                V_lower = V_ref(1:n);
                V_upper = V_ref(n+1:2*n);
                for i = 1:n
                    vol_range    = sys(i).max_V - sys(i).min_V;
                    track_err(i) = max([V_lower(i) - V_out(i), ...
                                        V_out(i)   - V_upper(i), 0]) / vol_range;
                end
        end

        fprintf('[t=%3d]  Power: %7.3f  Spill: %6.3f  TrackErr: %.4f\n', ...
            t, sum(p_out), sum(sp_out), mean(track_err));

        obj = sum(p_out);

        for i = 1:n
            X(i,:) = [V_out(i), p_out(i), u_out(i), sp_out(i), q_t(i)];
        end

    else
        warning('[t=%d] Gurobi status: %s. Returning NaN obj and zero X.', ...
            t, result.status);
    end
end

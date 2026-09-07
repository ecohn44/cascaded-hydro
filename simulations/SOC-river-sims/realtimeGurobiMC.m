function [result, obj, X, std_hat] = realtimeGurobiMC(t, c, eps, I_prev, q_error, V_prev, u_prev, V_ref, theta, lag, up_release, sys, model, bounds, framework, tracking)
% =========================================================================
% realtimeGurobiMC  — McCormick relaxation of the bilinear head term
%
% Identical to realtimeGurobi EXCEPT for C4 (power production constraint).
%
% CHANGE SUMMARY vs. realtimeGurobi
% ----------------------------------
%   1. nVars : 5*n  -->  6*n   (one extra z(i) per unit, block 5)
%   2. idx_z : new index function for z(i) = V(i)^b
%   3. buildVariables  : adds bounds block B6 for z(i)
%   4. buildLinearConstraints : replaces the fixed-head equality C4 with
%      four McCormick inequalities (MC1–MC4) using z(i) and u(i)
%   5. realtimeGurobiMC  : builds genconpow linking V(i) -> z(i) = V(i)^b
%   6. solveModel : accepts genconpow, sets FuncNonlinear = 1
%
% HEAD FORMULA (consistent with oracleGurobi)
%   h_i = a_i * V_i^{b_i}      (requires sys(i).a field)
%   z_i = V_i^{b_i}            (auxiliary, enforced by genconpow)
%   p_i = c * a_i * z_i * u_i  (bilinear, relaxed by McCormick)
%
% McCormick envelope of w = z*u over [zL,zU] x [uL,uU]:
%   MC1:  p/(c*a) >= zL*u + uL*z - zL*uL
%   MC2:  p/(c*a) >= zU*u + uU*z - zU*uU
%   MC3:  p/(c*a) <= zU*u + uL*z - zU*uL
%   MC4:  p/(c*a) <= zL*u + uU*z - zL*uU
%
% NOTE: The McCormick envelope is a RELAXATION — p may be overestimated
%       at non-corner solutions. This file is intended for testing only.
%       For production use, see realtimeGurobiPWL.m.
%
% INPUTS  (unchanged from realtimeGurobi)
%   t, c, eps, I_prev, q_error, V_prev, u_prev, V_ref, theta, lag,
%   up_release, sys, model, bounds, framework, tracking
%
% ADDED sys FIELD REQUIREMENT
%   sys(i).a  : head coefficient  (h = a * V^b)
%               if not present, computed as  a = min_h / min_V^b
% =========================================================================

    n = numel(sys);

    % --- derive sys(i).a if not provided ----------------------------
    for i = 1:n
        if ~isfield(sys(i), 'a') || isempty(sys(i).a)
            sys(i).a = sys(i).min_h / (sys(i).min_V ^ sys(i).b);
        end
    end

    % 1: Forecast inflow and estimate error (UNCHANGED)
    q_t = I_prev(:);
    if t > lag
        for i = 2:n
            q_t(i) = model.coef0 + model.coef1 * I_prev(i) + model.coef2 * up_release(i-1);
        end
    end
    std_hat = forecast_error(t, q_error, up_release, framework, model, sys);

    % 2: Calculate volume shift (UNCHANGED)
    switch bounds
        case {"det", "jcc-ssh"}
            for i = 1:n
                sys(i).V_eff_max = sys(i).max_V;
                sys(i).V_eff_min = sys(i).min_V;
            end
        case {"jcc-bon"}
            z_score = norminv(1 - (eps / (2*n)));
            V_min_shift =  z_score .* std_hat;
            V_max_shift = -z_score .* std_hat;
            for i = 1:n
                sys(i).V_eff_max = sys(i).max_V + V_max_shift(i);
                sys(i).V_eff_min = max(0, sys(i).min_V + V_min_shift(i)); % guard >= 0 for genconpow
            end
    end

    % 3: Decision variable layout
    % CHANGED: nVars 5*n -> 6*n, added idx_z
    nVars  = 6 * n;
    idx_V  = @(i)  (i - 1);
    idx_p  = @(i)  n   + (i - 1);
    idx_u  = @(i)  2*n + (i - 1);
    idx_sp = @(i)  3*n + (i - 1);
    idx_d  = @(i)  4*n + (i - 1);
    idx_z  = @(i)  5*n + (i - 1);   % NEW: z(i) = V(i)^b

    % 4: Variable bounds and objective coefficients
    [lb, ub, obj_coeff] = buildVariables(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_z, theta, sys);

    % 5: Linear constraints (C4 replaced with McCormick)
    [A, rhs, sense] = buildLinearConstraints(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_z, V_prev, u_prev, q_t, V_ref, c, sys);

    % 6: Build genconpow: z(i) = V(i)^{b_i}
    % This is what links the V decision variable to z via the nonlinear map.
    % Gurobi enforces this exactly via spatial branch-and-bound.
    genconpow = struct('xvar', {}, 'yvar', {}, 'a', {});
    for i = 1:n
        genconpow(i).xvar = idx_V(i) + 1;   % base  : V(i)
        genconpow(i).yvar = idx_z(i) + 1;   % result: z(i) = V(i)^a
        genconpow(i).a    = sys(i).b;        % exponent b_i in (0,1)
    end

    % 7: Solve
    result = solveModel(nVars, obj_coeff, A, rhs, sense, lb, ub, genconpow, t);

    % 8: Extract solution (UNCHANGED logic, no z reported externally)
    [obj, X] = extractSolution(result, t, n, idx_V, idx_p, idx_u, idx_sp, q_t, V_ref, sys, tracking);

end


% =========================================================================
% Helper: variable bounds and objective coefficients
% CHANGED: added idx_z argument, added block B6
% =========================================================================
function [lb, ub, obj_coeff] = buildVariables(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_z, theta, sys)

    lb        = zeros(nVars, 1);
    ub        =  inf(nVars, 1);
    obj_coeff =  zeros(nVars, 1);

    for i = 1:n

        % B1: Storage V(i)  — UNCHANGED
        lb(idx_V(i)+1)  = sys(i).V_eff_min;
        ub(idx_V(i)+1)  = sys(i).V_eff_max;

        % B2: Power output p(i)  — UNCHANGED
        lb(idx_p(i)+1)  = 0;
        ub(idx_p(i)+1)  = sys(i).F;
        obj_coeff(idx_p(i)+1) = -1;

        % B3: Turbine release u(i)  — UNCHANGED
        if t == 1
            lb(idx_u(i)+1) = sys(i).min_ut;
            ub(idx_u(i)+1) = sys(i).min_ut;
        else
            lb(idx_u(i)+1) = sys(i).min_ut;
            ub(idx_u(i)+1) = sys(i).max_ut;
        end

        % B4: Spill sp(i)  — UNCHANGED
        lb(idx_sp(i)+1) = 0;
        ub(idx_sp(i)+1) = inf;
        obj_coeff(idx_sp(i)+1) = 1e-4;

        % B5: Tracking error d(i)  — UNCHANGED
        lb(idx_d(i)+1)  = 0;
        ub(idx_d(i)+1)  = inf;
        obj_coeff(idx_d(i)+1) = theta / (sys(i).max_V - sys(i).min_V);

        % B6: NEW — auxiliary z(i) = V(i)^b
        % Since b in (0,1) and the power function is monotone increasing,
        % bounds map directly from the V bounds.
        % V_eff_min >= 0 is enforced in the bounds-shift block above.
        lb(idx_z(i)+1)  = sys(i).V_eff_min ^ sys(i).b;
        ub(idx_z(i)+1)  = sys(i).V_eff_max ^ sys(i).b;

    end
end


% =========================================================================
% Helper: linear constraints
% CHANGED: C4 replaced by four McCormick inequalities MC1–MC4
%          All other constraints (C1–C3, C5–C6) are IDENTICAL
% =========================================================================
function [A, rhs, sense] = buildLinearConstraints(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_z, V_prev, u_prev, q_t, V_ref, c, sys)

    col = @(x) x(:);   % force column vector — prevents vertcat errors

    rows  = zeros(0,1);
    cols  = zeros(0,1);
    vals  = zeros(0,1);
    rhs   = zeros(0,1);
    sense = '';
    row   = 0;

    for i = 1:n

        % (C1) Mass balance: V(i) + u(i) + sp(i) = V_prev(i) + q_t(i)
        % UNCHANGED
        row = row + 1;
        rows  = [rows;  row;          row;          row          ];
        cols  = [cols;  idx_V(i)+1;   idx_u(i)+1;   idx_sp(i)+1  ];
        vals  = [vals;  1;            1;            1            ];
        rhs   = [rhs;   V_prev(i) + q_t(i)                       ];
        sense = [sense; '='                                      ];

        if t > 1
            % (C2) Ramp-rate lower bound  — UNCHANGED
            row = row + 1;
            rows  = [rows;  row        ];
            cols  = [cols;  idx_u(i)+1 ];
            vals  = [vals;  -1         ];
            rhs   = [rhs;   -(u_prev(i) + sys(i).RR_dn)];
            sense = [sense; '<'        ];

            % (C3) Ramp-rate upper bound  — UNCHANGED
            row = row + 1;
            rows  = [rows;  row        ];
            cols  = [cols;  idx_u(i)+1 ];
            vals  = [vals;  1          ];
            rhs   = [rhs;   u_prev(i) + sys(i).RR_up];
            sense = [sense; '<'        ];
        end

        % ------------------------------------------------------------------
        % (C4) Power production — CHANGED
        %
        % Original (fixed head, pure LP):
        %   p(i) = c * h_fixed * u(i)   where h_fixed = h(V_prev)
        %
        % Now (McCormick relaxation using CURRENT V):
        %   p(i) = c * a_i * z(i) * u(i)
        %   z(i) = V(i)^{b_i}  (enforced by genconpow, not here)
        %
        % Bilinear product  w = z * u  over  z in [zL, zU], u in [uL, uU]
        % McCormick envelope (four linear inequalities):
        %   MC1:  w >= zL*u + uL*z - zL*uL   =>  -w + zL*u + uL*z <= zL*uL
        %   MC2:  w >= zU*u + uU*z - zU*uU   =>  -w + zU*u + uU*z <= zU*uU
        %   MC3:  w <= zU*u + uL*z - zU*uL   =>   w - zU*u - uL*z <= -zU*uL
        %   MC4:  w <= zL*u + uU*z - zL*uU   =>   w - zL*u - uU*z <= -zL*uU
        %
        % Substituting  w = p / (c * a_i):
        % ------------------------------------------------------------------
        zL  = sys(i).V_eff_min ^ sys(i).b;   % lower bound on z(i)
        zU  = sys(i).V_eff_max ^ sys(i).b;   % upper bound on z(i)
        uL  = sys(i).min_ut;                  % lower bound on u(i)
        uU  = sys(i).max_ut;                  % upper bound on u(i)
        ca  = c * sys(i).a;                   % combined scaling factor

        % MC1: -p/(ca) + zL*u + uL*z <= zL*uL
        row = row + 1;
        rows  = [rows;  row;          row;          row          ];
        cols  = [cols;  idx_p(i)+1;   idx_u(i)+1;   idx_z(i)+1   ];
        vals  = [vals;  -1/ca;        zL;           uL           ];
        rhs   = [rhs;   zL * uL                                   ];
        sense = [sense; '<'                                       ];

        % MC2: -p/(ca) + zU*u + uU*z <= zU*uU
        row = row + 1;
        rows  = [rows;  row;          row;          row          ];
        cols  = [cols;  idx_p(i)+1;   idx_u(i)+1;   idx_z(i)+1   ];
        vals  = [vals;  -1/ca;        zU;           uU           ];
        rhs   = [rhs;   zU * uU                                   ];
        sense = [sense; '<'                                       ];

        % MC3: p/(ca) - zU*u - uL*z <= -zU*uL
        row = row + 1;
        rows  = [rows;  row;          row;          row          ];
        cols  = [cols;  idx_p(i)+1;   idx_u(i)+1;   idx_z(i)+1   ];
        vals  = [vals;  1/ca;        -zU;          -uL           ];
        rhs   = [rhs;  -zU * uL                                   ];
        sense = [sense; '<'                                       ];

        % MC4: p/(ca) - zL*u - uU*z <= -zL*uU
        row = row + 1;
        rows  = [rows;  row;          row;          row          ];
        cols  = [cols;  idx_p(i)+1;   idx_u(i)+1;   idx_z(i)+1   ];
        vals  = [vals;  1/ca;        -zL;          -uU           ];
        rhs   = [rhs;  -zL * uU                                   ];
        sense = [sense; '<'                                       ];

        % (C5) Tracking error upper  — UNCHANGED
        V_upper = V_ref(n+1:2*n);
        row = row + 1;
        rows  = [rows;  row;          row          ];
        cols  = [cols;  idx_V(i)+1;   idx_d(i)+1   ];
        vals  = [vals;  1;           -1            ];
        rhs   = [rhs;   V_upper(i)                  ];
        sense = [sense; '<'                         ];

        % (C6) Tracking error lower  — UNCHANGED
        V_lower = V_ref(1:n);
        row = row + 1;
        rows  = [rows;  row;          row          ];
        cols  = [cols;  idx_V(i)+1;   idx_d(i)+1   ];
        vals  = [vals;  -1;          -1            ];
        rhs   = [rhs;  -V_lower(i)                  ];
        sense = [sense; '<'                         ];

    end

    nCons = row;
    A = sparse(col(rows), col(cols), col(vals), nCons, nVars);
end


% =========================================================================
% Helper: solve
% CHANGED: added genconpow argument, FuncNonlinear = 1
% =========================================================================
function result = solveModel(nVars, obj_coeff, A, rhs, sense, lb, ub, genconpow, t_period)

    grb_model.modelname  = sprintf('hydroRT_MC_t%d', t_period);
    grb_model.modelsense = 'min';
    grb_model.obj        = obj_coeff(:);
    grb_model.A          = A;
    grb_model.rhs        = rhs(:);
    grb_model.sense      = sense;
    grb_model.lb         = lb(:);
    grb_model.ub         = ub(:);
    grb_model.vtype      = repmat('C', nVars, 1);

    % Attach genconpow: z(i) = V(i)^{b_i}
    if ~isempty(genconpow)
        grb_model.genconpow = genconpow;
    end

    params.OutputFlag    = 0;
    params.Seed          = 1;
    params.TimeLimit     = 10;
    params.Threads       = 1;    % fixed for deterministic real-time latency
    params.FuncNonlinear = 1;    % NEW: activates spatial B&B for genconpow

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
% Helper: extract solution  — UNCHANGED from realtimeGurobi
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

        % Tracking error  — UNCHANGED
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

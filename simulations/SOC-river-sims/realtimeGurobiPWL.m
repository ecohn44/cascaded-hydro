function [result, obj, X, std_hat] = realtimeGurobiPWL(t, c, eps, I_prev, q_error, V_prev, u_prev, V_ref, theta, lag, up_release, sys, model, bounds, framework, tracking)
% =========================================================================
% realtimeGurobiPWL  —  Real-time single-period hydropower dispatch
%                       Piecewise-Linear MILP formulation
%
% POWER PRODUCTION MODEL
% ----------------------
%   p(i) = c * h_bar(k*) * u(i)
%
%   where k* is the segment containing the CURRENT V(i), and
%   h_bar(k) = a_i * V_mid_k ^ b_i  (precomputed scalar per segment).
%
%   This is a pure MILP — no bilinear terms, no McCormick relaxation gap.
%
%   1. Binary y(i,k) in {0,1} selects the segment k containing V(i).
%   2. v(i,k) = y(i,k) * u(i) is EXACTLY linearised by four Big-M
%      constraints.  Because y is binary (not continuous), there is
%      no relaxation gap at any integer-feasible node.
%   3. p(i) = c * SUM_k h_bar(k) * v(i,k)  — fully linear.
%
% VARIABLE LAYOUT  (5*n + 2*n*K variables total)
% -----------------------------------------------
%   Block 0  [0 .. n-1]            V(i)      storage
%   Block 1  [n .. 2n-1]           p(i)      power
%   Block 2  [2n .. 3n-1]          u(i)      release
%   Block 3  [3n .. 4n-1]          sp(i)     spill
%   Block 4  [4n .. 5n-1]          d(i)      tracking slack
%   Block 5  [5n .. 5n+nK-1]       y(i,k)    binary segment
%   Block 6  [5n+nK .. 5n+2nK-1]   v(i,k)    auxiliary  v = y*u
%
% ACCURACY CONTROL
% ----------------
%   K  (segments per unit) trades approximation accuracy for MILP size.
%   Default K = 10.  Override per-unit via sys(i).K before calling.
%
%   Approximation error bound per segment:
%     eps_K <= (c * u_max_i / 2) * max_k |V_bp(k+1)^b - V_bp(k)^b|
%
% INPUTS  (interface unchanged from original realtimeGurobi)
%   t          : Current time period index
%   c          : Power conversion coefficient (scalar)
%   eps        : JCC violation probability
%   I_prev     : Inflow vector (n x 1)
%   q_error    : Forecast error history struct
%   V_prev     : Measured storage at start of period (n x 1)
%   u_prev     : Turbine release in previous period (n x 1)
%   V_ref      : [V_lower(1..n); V_upper(1..n)]  (2n x 1)
%   theta      : Volume-tracking penalty weight (scalar)
%   lag        : Warm-up lag before regression forecast activates
%   up_release : Upstream release vector
%   sys        : Struct array — fields: min_V, max_V, min_h, max_h,
%                  a, b, min_ut, max_ut, RR_dn, RR_up, F, K (optional)
%   model      : Inflow forecast regression model struct
%   bounds     : 'det' | 'jcc-ssh' | 'jcc-bon'
%   framework  : Forecast framework identifier
%   tracking   : 'mean' | 'envelope'
%
% OUTPUTS
%   result   : Gurobi result struct
%   obj      : Total PHYSICAL power generation (scalar)
%   X        : [V, p_physical, u, sp, q] per unit  (n x 5)
%   std_hat  : Forecast standard deviation  (n x 1)
% =========================================================================

    n = numel(sys);

    % Default number of PWL segments — overridable via sys(i).K
    K = 10;
    if isfield(sys(1), 'K') && ~isempty(sys(1).K)
        K = sys(1).K;
    end

    % ------------------------------------------------------------------
    % 1: Forecast inflow
    % ------------------------------------------------------------------
    q_t = I_prev(:);
    if t > lag
        for i = 2:n
            q_t(i) = model.coef0 + model.coef1 * I_prev(i) ...
                                  + model.coef2 * up_release(i-1);
        end
    end
    std_hat = forecast_error(t, q_error, up_release, framework, model, sys);

    % ------------------------------------------------------------------
    % 2: Effective volume bounds
    % ------------------------------------------------------------------
    switch bounds
        case {"det", "jcc-ssh"}
            for i = 1:n
                sys(i).V_eff_max = sys(i).max_V;
                sys(i).V_eff_min = sys(i).min_V;
            end

        case {"jcc-bon"}
            z_score     = norminv(1 - (eps / (2*n)));
            V_min_shift =  z_score .* std_hat;
            V_max_shift = -z_score .* std_hat;
            for i = 1:n
                % Guard: genconpow and V^b require V_eff_min >= 0
                sys(i).V_eff_min = max(0, sys(i).min_V + V_min_shift(i));
                sys(i).V_eff_max = sys(i).max_V + V_max_shift(i);
            end
    end

    % ------------------------------------------------------------------
    % 3: Index functions  (all return 0-based; add 1 for MATLAB indexing)
    % ------------------------------------------------------------------
    nVars  = 5*n + 2*n*K;

    idx_V  = @(i)      (i - 1);
    idx_p  = @(i)      n           + (i - 1);
    idx_u  = @(i)      2*n         + (i - 1);
    idx_sp = @(i)      3*n         + (i - 1);
    idx_d  = @(i)      4*n         + (i - 1);
    idx_y  = @(i,k)    5*n         + (i-1)*K + (k-1);
    idx_v  = @(i,k)    5*n + n*K   + (i-1)*K + (k-1);

    % ------------------------------------------------------------------
    % 4: Variable bounds, types and objective coefficients
    % ------------------------------------------------------------------
    [lb, ub, obj_coeff, vtype] = buildVariables( ...
        t, n, K, nVars, ...
        idx_V, idx_p, idx_u, idx_sp, idx_d, idx_y, idx_v, ...
        theta, sys);

    % ------------------------------------------------------------------
    % 5: Linear constraint matrix
    % ------------------------------------------------------------------
    [A, rhs, sense] = buildLinearConstraints( ...
        t, n, K, nVars, ...
        idx_V, idx_p, idx_u, idx_sp, idx_d, idx_y, idx_v, ...
        V_prev, u_prev, q_t, V_ref, c, sys);

    % ------------------------------------------------------------------
    % 6: Solve
    % ------------------------------------------------------------------
    result = solveModel(nVars, obj_coeff, A, rhs, sense, lb, ub, vtype, t);

    % ------------------------------------------------------------------
    % 7: Extract solution
    % ------------------------------------------------------------------
    [obj, X] = extractSolution( ...
        result, t, n, K, ...
        idx_V, idx_p, idx_u, idx_sp, idx_y, idx_v, ...
        q_t, V_ref, c, sys, tracking);

end


% =========================================================================
%  HELPER: variable bounds, types and objective coefficients
% =========================================================================
function [lb, ub, obj_coeff, vtype] = buildVariables( ...
        t, n, K, nVars, ...
        idx_V, idx_p, idx_u, idx_sp, idx_d, idx_y, idx_v, ...
        theta, sys)

    lb        = zeros(nVars, 1);
    ub        = inf(nVars, 1);
    obj_coeff = zeros(nVars, 1);
    vtype     = repmat('C', nVars, 1);

    for i = 1:n

        % B1: Storage
        lb(idx_V(i)+1) = sys(i).V_eff_min;
        ub(idx_V(i)+1) = sys(i).V_eff_max;

        % B2: Power
        lb(idx_p(i)+1)        = 0;
        ub(idx_p(i)+1)        = sys(i).F;
        obj_coeff(idx_p(i)+1) = -1;

        % B3: Turbine release  (pin at min for t=1 initial condition)
        if t == 1
            lb(idx_u(i)+1) = sys(i).min_ut;
            ub(idx_u(i)+1) = sys(i).min_ut;
        else
            lb(idx_u(i)+1) = sys(i).min_ut;
            ub(idx_u(i)+1) = sys(i).max_ut;
        end

        % B4: Spill
        lb(idx_sp(i)+1)        = 0;
        ub(idx_sp(i)+1)        = inf;
        obj_coeff(idx_sp(i)+1) = 1e-4;

        % B5: Tracking error slack
        lb(idx_d(i)+1)        = 0;
        ub(idx_d(i)+1)        = inf;
        vol_range              = sys(i).max_V - sys(i).min_V;
        obj_coeff(idx_d(i)+1) = theta / vol_range;

        % B6: Binary segment indicators y(i,k)
        for k = 1:K
            lb(idx_y(i,k)+1)    = 0;
            ub(idx_y(i,k)+1)    = 1;
            vtype(idx_y(i,k)+1) = 'B';
        end

        % B7: Auxiliary continuous  v(i,k) = y(i,k)*u(i)
        for k = 1:K
            lb(idx_v(i,k)+1) = 0;
            ub(idx_v(i,k)+1) = sys(i).max_ut;
        end

    end
end


% =========================================================================
%  HELPER: sparse linear constraint matrix
%
%  ORIENTATION CONVENTION
%  ----------------------
%  rows, cols, vals are ALWAYS (N x 1) column vectors.
%
%  The local helper   col = @(x) x(:)   enforces this at every append,
%  making orientation errors structurally impossible regardless of how
%  MATLAB evaluates intermediate expressions (scalar, row, or column).
%
%  Rule: wrap EVERY non-scalar expression appended to rows/cols/vals
%  with col(...) or col(scalar) — scalars are coerced to 1x1 by (:).
% =========================================================================
function [A, rhs, sense] = buildLinearConstraints( ...
        t, n, K, nVars, ...
        idx_V, idx_p, idx_u, idx_sp, idx_d, idx_y, idx_v, ...
        V_prev, u_prev, q_t, V_ref, c, sys)

    % ---------------------------------------------------------------
    % Orientation helper: guarantees any expression becomes (M x 1).
    % col(scalar)  -> [scalar]  1x1
    % col(row_vec) -> column    Kx1
    % col(col_vec) -> unchanged Kx1
    % This is the ONLY defensive measure needed — apply consistently.
    % ---------------------------------------------------------------
    col = @(x) x(:);

    rows  = zeros(0,1);
    cols  = zeros(0,1);
    vals  = zeros(0,1);
    rhs   = zeros(0,1);
    sense = '';
    row   = 0;

    for i = 1:n

        % ----------------------------------------------------------
        % Segment geometry for unit i
        % ----------------------------------------------------------
        V_lo = sys(i).V_eff_min;
        V_hi = sys(i).V_eff_max;
        dV   = (V_hi - V_lo) / K;

        % Breakpoints  V_bp(k) = V_lo + (k-1)*dV,  k = 1..K+1
        % Declared explicitly as (K+1) x 1 column vector.
        V_bp = V_lo + col(0:K) * dV;            % (K+1) x 1

        % Midpoint head values  h_bar(k) = a_i * V_mid_k ^ b_i
        % Declared explicitly as K x 1 column vector.
        h_bar = zeros(K, 1);
        for k = 1:K
            V_mid_k  = 0.5 * (V_bp(k) + V_bp(k+1));
            h_bar(k) = sys(i).a * V_mid_k ^ sys(i).b;
        end
        % Defensive coercion: guarantee column even if K==1
        h_bar = col(h_bar);                      % K x 1

        u_min_i = sys(i).min_ut;
        u_max_i = sys(i).max_ut;

        % ==========================================================
        % (C1) Mass balance
        %      V(i) + u(i) + sp(i) = V_prev(i) + q_t(i)
        % ==========================================================
        row   = row + 1;
        rows  = [rows;  col([row; row; row])                                  ];
        cols  = [cols;  col([idx_V(i)+1; idx_u(i)+1; idx_sp(i)+1])           ];
        vals  = [vals;  col([1; 1; 1])                                        ];
        rhs   = [rhs;   V_prev(i) + q_t(i)                                   ];
        sense = [sense; '='                                                   ];

        % ==========================================================
        % (C2) Ramp-rate lower: u(i) >= u_prev(i) - RR_dn
        %      written as  -u(i) <= -(u_prev(i) + RR_dn)  [note sign]
        % (C3) Ramp-rate upper: u(i) <= u_prev(i) + RR_up
        % ==========================================================
        if t > 1
            row   = row + 1;
            rows  = [rows;  row                         ];
            cols  = [cols;  idx_u(i)+1                  ];
            vals  = [vals;  -1                          ];
            rhs   = [rhs;   -(u_prev(i) + sys(i).RR_dn)];
            sense = [sense; '<'                         ];

            row   = row + 1;
            rows  = [rows;  row                        ];
            cols  = [cols;  idx_u(i)+1                 ];
            vals  = [vals;  1                          ];
            rhs   = [rhs;   u_prev(i) + sys(i).RR_up  ];
            sense = [sense; '<'                        ];
        end

        % ==========================================================
        % (C4) PWL-MILP power production block
        %
        % GOAL:  p(i) = c * a_i * V(i)^b_i * u(i)
        %
        % Step 1 — Segment selection [C4a, C4b, C4c]
        %   y(i,k) in {0,1},  SUM_k y(i,k) = 1         [C4a]
        %   y(i,k)=1 => V(i) in [V_bp(k), V_bp(k+1)]   [C4b, C4c]
        %
        % Step 2 — Head approximation
        %   h(i) ~ h_bar(k*)  where h_bar(k) is a precomputed scalar.
        %
        % Step 3 — Exact binary linearisation of v(i,k) = y(i,k)*u(i)
        %   Because y is BINARY:
        %     y=0 => C4d + C4g force v(i,k) = 0  exactly
        %     y=1 => C4e + C4f force v(i,k) = u(i) exactly
        %   This is NOT an approximation.
        %
        % Step 4 — Power and release equalities [C4h, C4i]
        %   p(i) = c * SUM_k h_bar(k) * v(i,k)
        %   u(i) = SUM_k v(i,k)          [tightens LP relaxation]
        % ==========================================================

        % ----------------------------------------------------------
        % C4a: SUM_k y(i,k) = 1
        %
        % y_cols and repmat(row, K, 1) are explicitly K x 1 columns.
        % col() wraps the arrayfun output for safety when K==1.
        % ----------------------------------------------------------
        row    = row + 1;
        y_cols = col(arrayfun(@(k) idx_y(i,k)+1, 1:K));   % K x 1
        rows   = [rows;  col(repmat(row, K, 1))           ];
        cols   = [cols;  y_cols                            ];
        vals   = [vals;  col(ones(K, 1))                   ];
        rhs    = [rhs;   1                                 ];
        sense  = [sense; '='                               ];

        % ----------------------------------------------------------
        % Per-segment constraints C4b–C4g
        % ----------------------------------------------------------
        for k = 1:K

            % Tight Big-M values derived from segment geometry.
            % Using problem bounds (not a large constant) keeps the
            % LP relaxation as tight as possible at every B&B node.
            M_lo = (k-1) * dV;    % distance V_lo -> V_bp(k)
            M_hi = (K-k) * dV;    % distance V_bp(k+1) -> V_hi

            % C4b: -V(i) + M_lo * y(i,k) <= -V_lo
            %      Enforces V(i) >= V_bp(k) when y(i,k) = 1
            row   = row + 1;
            rows  = [rows;  col([row; row])                         ];
            cols  = [cols;  col([idx_V(i)+1; idx_y(i,k)+1])        ];
            vals  = [vals;  col([-1; M_lo])                         ];
            rhs   = [rhs;   -V_lo                                   ];
            sense = [sense; '<'                                     ];

            % C4c:  V(i) + M_hi * y(i,k) <= V_hi
            %      Enforces V(i) <= V_bp(k+1) when y(i,k) = 1
            row   = row + 1;
            rows  = [rows;  col([row; row])                         ];
            cols  = [cols;  col([idx_V(i)+1; idx_y(i,k)+1])        ];
            vals  = [vals;  col([1; M_hi])                          ];
            rhs   = [rhs;   V_hi                                    ];
            sense = [sense; '<'                                     ];

            % C4d: v(i,k) >= u_min * y(i,k)
            %      -v(i,k) + u_min * y(i,k) <= 0
            row   = row + 1;
            rows  = [rows;  col([row; row])                         ];
            cols  = [cols;  col([idx_v(i,k)+1; idx_y(i,k)+1])      ];
            vals  = [vals;  col([-1; u_min_i])                      ];
            rhs   = [rhs;   0                                       ];
            sense = [sense; '<'                                     ];

            % C4e: v(i,k) >= u(i) - u_max*(1-y(i,k))
            %      -v(i,k) + u(i) + u_max*y(i,k) <= u_max
            row   = row + 1;
            rows  = [rows;  col([row; row; row])                               ];
            cols  = [cols;  col([idx_v(i,k)+1; idx_u(i)+1; idx_y(i,k)+1])     ];
            vals  = [vals;  col([-1; 1; u_max_i])                              ];
            rhs   = [rhs;   u_max_i                                            ];
            sense = [sense; '<'                                                ];

            % C4f: v(i,k) <= u(i) - u_min*(1-y(i,k))
            %       v(i,k) - u(i) - u_min*y(i,k) <= -u_min
            row   = row + 1;
            rows  = [rows;  col([row; row; row])                               ];
            cols  = [cols;  col([idx_v(i,k)+1; idx_u(i)+1; idx_y(i,k)+1])     ];
            vals  = [vals;  col([1; -1; -u_min_i])                             ];
            rhs   = [rhs;   -u_min_i                                           ];
            sense = [sense; '<'                                                ];

            % C4g: v(i,k) <= u_max * y(i,k)
            %       v(i,k) - u_max * y(i,k) <= 0
            row   = row + 1;
            rows  = [rows;  col([row; row])                         ];
            cols  = [cols;  col([idx_v(i,k)+1; idx_y(i,k)+1])      ];
            vals  = [vals;  col([1; -u_max_i])                      ];
            rhs   = [rhs;   0                                       ];
            sense = [sense; '<'                                     ];

        end  % for k = 1:K

        % ----------------------------------------------------------
        % C4h: Power equality
        %      p(i) - c * SUM_k h_bar(k)*v(i,k) = 0
        %
        % v_cols is K x 1 (col applied).
        % -c * h_bar is K x 1 (h_bar already coerced above).
        % col() applied to both for guaranteed orientation.
        % ----------------------------------------------------------
        v_cols = col(arrayfun(@(k) idx_v(i,k)+1, 1:K));   % K x 1
        row    = row + 1;
        rows   = [rows;  row;                col(repmat(row, K, 1)) ];
        cols   = [cols;  idx_p(i)+1;         v_cols                 ];
        vals   = [vals;  1;                  col(-c .* h_bar)       ];
        rhs    = [rhs;   0                                          ];
        sense  = [sense; '='                                        ];

        % ----------------------------------------------------------
        % C4i: Release consistency
        %      u(i) - SUM_k v(i,k) = 0
        %      Redundant at integer nodes but tightens LP relaxation.
        % ----------------------------------------------------------
        row    = row + 1;
        rows   = [rows;  row;                col(repmat(row, K, 1)) ];
        cols   = [cols;  idx_u(i)+1;         v_cols                 ];
        vals   = [vals;  1;                  col(-ones(K, 1))       ];
        rhs    = [rhs;   0                                          ];
        sense  = [sense; '='                                        ];

        % ==========================================================
        % (C5) Tracking upper: V(i) - d(i) <= V_upper_ref(i)
        % (C6) Tracking lower: -V(i) - d(i) <= -V_lower_ref(i)
        % ==========================================================
        V_upper = V_ref(n+1:2*n);
        V_lower = V_ref(1:n);

        row   = row + 1;
        rows  = [rows;  col([row; row])                     ];
        cols  = [cols;  col([idx_V(i)+1; idx_d(i)+1])       ];
        vals  = [vals;  col([1; -1])                         ];
        rhs   = [rhs;   V_upper(i)                          ];
        sense = [sense; '<'                                 ];

        row   = row + 1;
        rows  = [rows;  col([row; row])                     ];
        cols  = [cols;  col([idx_V(i)+1; idx_d(i)+1])       ];
        vals  = [vals;  col([-1; -1])                        ];
        rhs   = [rhs;   -V_lower(i)                         ];
        sense = [sense; '<'                                 ];

    end  % for i = 1:n

    nCons = row;
    A = sparse(rows, cols, vals, nCons, nVars);
end


% =========================================================================
%  HELPER: configure and solve the Gurobi MILP
% =========================================================================
function result = solveModel(nVars, obj_coeff, A, rhs, sense, lb, ub, vtype, t_period)

    grb_model.modelname  = sprintf('hydroRT_pwl_t%d', t_period);
    grb_model.modelsense = 'min';
    grb_model.obj        = obj_coeff(:);    % enforce column
    grb_model.A          = A;
    grb_model.rhs        = rhs(:);          % enforce column
    grb_model.sense      = sense;
    grb_model.lb         = lb(:);           % enforce column
    grb_model.ub         = ub(:);           % enforce column
    grb_model.vtype      = vtype(:);        % enforce column
    % Pure MILP: no genconpow, no FuncNonlinear required.

    params.OutputFlag = 0;
    params.Seed       = 1;
    params.TimeLimit  = 10;
    params.Threads    = 1;
    params.MIPGap     = 1e-4;

    result = gurobi(grb_model, params);

    if strcmp(result.status, 'INFEASIBLE')
        fprintf('[t=%d] INFEASIBLE — running IIS diagnosis.\n', t_period);
        iis      = gurobi_iis(grb_model);
        iis_rows = find(iis.Arows);
        for r = iis_rows'
            fprintf('  Row %d | sense=%s | rhs=%.6f\n', ...
                r, grb_model.sense(r), grb_model.rhs(r));
            [~, cs, vs] = find(grb_model.A(r,:));
            for j = 1:numel(cs)
                fprintf('    col %d  coeff=%.4f  lb=%.4f  ub=%.4f\n', ...
                    cs(j), vs(j), grb_model.lb(cs(j)), grb_model.ub(cs(j)));
            end
        end
    end
end


% =========================================================================
%  HELPER: extract and report the solution
% =========================================================================
function [obj, X] = extractSolution( ...
        result, t, n, K, ...
        idx_V, idx_p, idx_u, idx_sp, idx_y, idx_v, ...
        q_t, V_ref, c, sys, tracking)

    obj = NaN;
    X   = zeros(n, 5);   % [V, p_physical, u, sp, q]

    valid_statuses = {'OPTIMAL', 'SUBOPTIMAL', 'TIME_LIMIT'};
    has_solution   = ismember(result.status, valid_statuses) ...
                     && isfield(result, 'x') && ~isempty(result.x);

    if ~has_solution
        warning('[t=%d] Gurobi status: %s. Returning NaN obj and zero X.', ...
            t, result.status);
        return
    end

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

    % ------------------------------------------------------------------
    % Physical power using exact nonlinear head at optimised V
    %   p_out      = MILP PWL approximation (may have small bias)
    %   p_physical = true physical output   (returned in X and obj)
    % ------------------------------------------------------------------
    p_physical = zeros(n,1);
    for i = 1:n
        h_exact       = sys(i).a * V_out(i) ^ sys(i).b;
        p_physical(i) = c * h_exact * u_out(i);
    end

    % ------------------------------------------------------------------
    % Console diagnostics
    % ------------------------------------------------------------------
    p_pwl_err        = p_out - p_physical;
    [max_err, i_max] = max(abs(p_pwl_err));

    fprintf('[t=%3d]  Power(MILP): %7.3f  Power(exact): %7.3f  ', ...
        t, sum(p_out), sum(p_physical));
    fprintf('MaxPWLerr: %.4f (unit %d)  Spill: %.3f\n', ...
        max_err, i_max, sum(sp_out));

    for i = 1:n
        V_lo    = sys(i).V_eff_min;
        V_hi    = sys(i).V_eff_max;
        dV      = (V_hi - V_lo) / K;
        y_vals  = x(arrayfun(@(k) idx_y(i,k)+1, 1:K));
        active_k = find(y_vals > 0.5, 1);
        if ~isempty(active_k)
            V_mid_k = V_lo + (active_k - 0.5) * dV;
            h_pwl   = sys(i).a * V_mid_k ^ sys(i).b;
            h_exact = sys(i).a * V_out(i) ^ sys(i).b;
            fprintf('         unit %d: seg=%d/%d  h_pwl=%.4f  h_exact=%.4f  V=%.4f  u=%.4f\n', ...
                i, active_k, K, h_pwl, h_exact, V_out(i), u_out(i));
        end
    end

    % ------------------------------------------------------------------
    % Tracking error
    % ------------------------------------------------------------------
    track_err = zeros(n,1);
    switch tracking
        case "mean"
            V_mid_ref = (V_ref(1:n) + V_ref(n+1:2*n)) / 2;
            for i = 1:n
                vol_range    = sys(i).max_V - sys(i).min_V;
                track_err(i) = abs(V_out(i) - V_mid_ref(i)) / vol_range;
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

    fprintf('         TrackErr: %.4f\n', mean(track_err));

    % ------------------------------------------------------------------
    % Return physical values — not the MILP approximation
    % ------------------------------------------------------------------
    obj = sum(p_physical);

    for i = 1:n
        X(i,:) = [V_out(i), p_physical(i), u_out(i), sp_out(i), q_t(i)];
    end
end

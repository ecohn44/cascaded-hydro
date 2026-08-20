function [result, obj, X] = oracleGurobi(T, c, I, lag, s)
% ORACLE  Generalized n-unit Hydropower Optimization Model
%
% Solves the hydropower dispatch over T time steps using the Gurobi
% MATLAB API directly. Replaces YALMIP because YALMIP cannot translate
% the nonlinear head-storage relationship:
%
%   h(i,t) = a_i * V(i,t)^b_i   (b_i in (0,1), V(i,t) > 0)
%
% to Gurobi. Instead, we introduce an auxiliary variable:
%
%   z(i,t) = V(i,t)^b_i   [Gurobi POW general constraint]
%
% and substitute h(i,t) = a_i * z(i,t) into all linear constraints,
% eliminating h as a decision variable entirely.
%
% FuncNonlinear=1 activates Gurobi's spatial branch-and-bound, which
% provides an EXACT (globally optimal) solution.
%
% INPUTS
%   T   : Number of time periods in optimization horizon
%   c   : Power conversion coefficient
%   I   : Exogenous inflow matrix (n x T)
%   lag : Cascade travel lag (integer, time steps)
%   s   : Struct array (1 x n) of per-unit parameters:
%           s(i).a      - head coefficient
%           s(i).b      - head exponent (0 < b < 1)
%           s(i).V0     - initial storage volume
%           s(i).min_V  - minimum volume bound
%           s(i).max_V  - maximum volume bound
%           s(i).min_ut - minimum release
%           s(i).max_ut - maximum release
%           s(i).RR_dn  - ramp-down rate limit
%           s(i).RR_up  - ramp-up rate limit
%           s(i).F      - feeder (power) capacity
%
% OUTPUTS
%   result : Gurobi result struct
%   obj    : Total power generation over horizon (scalar)
%   X      : Solution matrix [V, p, u, sp, q] stacked per unit (T x 5n)

    n = numel(s);
    X = [];

    % =====================================================================
    %% Variable Index Helpers (0-based for Gurobi MATLAB API)
    %
    % Decision variables:
    %   V (i,t)  : reservoir storage volume         [n x T]
    %   p (i,t)  : power output                     [n x T]
    %   u (i,t)  : water release (turbine discharge) [n x T]
    %   sp(i,t)  : spill flow                        [n x T]
    %   z (i,t)  : auxiliary = V(i,t)^b_i, t=2..T  [n x (T-1)]
    %
    % Total variables: 4*n*T + n*(T-1)
    % =====================================================================
    nVars = 4*n*T + n*(T-1);

    idx_V  = @(i,t)  (i-1)*T         + (t-1);          % block 0
    idx_p  = @(i,t)  n*T              + (i-1)*T + (t-1); % block 1
    idx_u  = @(i,t)  2*n*T            + (i-1)*T + (t-1); % block 2
    idx_sp = @(i,t)  3*n*T            + (i-1)*T + (t-1); % block 3
    idx_z  = @(i,t)  4*n*T            + (i-1)*(T-1) + (t-2); % block 4, t=2..T only

    % =====================================================================
    %% Variable Bounds and Objective Coefficients
    % Gurobi minimizes: we minimize -sum(p) + 1e-4*sum(sp)
    % =====================================================================
    lb        = -inf(nVars, 1);
    ub        =  inf(nVars, 1);
    obj_coeff =  zeros(nVars, 1);

    for i = 1:n
        for t = 1:T
            % --- Volume: min_V <= V(i,t) <= max_V ---
            lb(idx_V(i,t)+1)  = s(i).min_V;
            ub(idx_V(i,t)+1)  = s(i).max_V;

            % --- Power: 0 <= p(i,t) <= F; objective: -p ---
            lb(idx_p(i,t)+1)  = 0;
            ub(idx_p(i,t)+1)  = s(i).F;
            obj_coeff(idx_p(i,t)+1) = -1;

            % --- Release: min_ut <= u(i,t) <= max_ut ---
            lb(idx_u(i,t)+1)  = s(i).min_ut;
            ub(idx_u(i,t)+1)  = s(i).max_ut;

            % --- Spill: 0 <= sp(i,t); objective: +1e-4*sp ---
            lb(idx_sp(i,t)+1) = 0;
            obj_coeff(idx_sp(i,t)+1) = 1e-4;

            % --- Auxiliary z(i,t) for t >= 2: bounds from V bounds ---
            if t >= 2
                lb(idx_z(i,t)+1) = s(i).min_V^s(i).b;
                ub(idx_z(i,t)+1) = s(i).max_V^s(i).b;
            end
        end
    end

    % =====================================================================
    %% Linear Constraints (assembled as sparse triplets)
    % =====================================================================
    rows  = [];
    cols  = [];
    vals  = [];
    rhs   = [];
    sense = '';
    row   = 0;

    for t = 1:T
        for i = 1:n

            if t == 1
                % ---------------------------------------------------------
                % Mass balance at t=1:
                %   V(i,1) + u(i,1) + sp(i,1) = V0 + I(i,1)
                % ---------------------------------------------------------
                row = row + 1;
                rows  = [rows;  row;            row;            row            ];
                cols  = [cols;  idx_V(i,1)+1;   idx_u(i,1)+1;  idx_sp(i,1)+1  ];
                vals  = [vals;  1;              1;              1              ];
                rhs   = [rhs;   s(i).V0 + I(i,1)];
                sense = [sense; '='];

                % ---------------------------------------------------------
                % Initial release fixed to minimum:
                %   u(i,1) = min_ut
                % ---------------------------------------------------------
                row = row + 1;
                rows  = [rows;  row           ];
                cols  = [cols;  idx_u(i,1)+1  ];
                vals  = [vals;  1             ];
                rhs   = [rhs;   s(i).min_ut   ];
                sense = [sense; '='];

                % ---------------------------------------------------------
                % Power at t=1 (head computed from constant V0):
                %   p(i,1) - c * h0 * u(i,1) = 0
                %   where h0 = a_i * V0^b_i  (scalar constant)
                % ---------------------------------------------------------
                h0 = s(i).a * (s(i).V0^s(i).b);
                row = row + 1;
                rows  = [rows;  row;            row           ];
                cols  = [cols;  idx_p(i,1)+1;   idx_u(i,1)+1  ];
                vals  = [vals;  1;              -c*h0         ];
                rhs   = [rhs;   0              ];
                sense = [sense; '='];

            else  % t >= 2

                % ---------------------------------------------------------
                % Ramp-rate lower bound: u(i,t) - u(i,t-1) >= RR_dn
                %   <==>  -u(i,t) + u(i,t-1) <= -RR_dn
                % ---------------------------------------------------------
                row = row + 1;
                rows  = [rows;  row;              row             ];
                cols  = [cols;  idx_u(i,t)+1;     idx_u(i,t-1)+1  ];
                vals  = [vals;  -1;               1               ];
                rhs   = [rhs;   -s(i).RR_dn       ];
                sense = [sense; '<'];

                % ---------------------------------------------------------
                % Ramp-rate upper bound: u(i,t) - u(i,t-1) <= RR_up
                % ---------------------------------------------------------
                row = row + 1;
                rows  = [rows;  row;              row             ];
                cols  = [cols;  idx_u(i,t)+1;     idx_u(i,t-1)+1  ];
                vals  = [vals;  1;                -1              ];
                rhs   = [rhs;   s(i).RR_up        ];
                sense = [sense; '<'];

                % ---------------------------------------------------------
                % McCormick Envelopes for p(i,t) = c * h(i,t) * u(i,t)
                %
                % Since h(i,t) = a_i * z(i,t), we substitute directly:
                %   h_min = a_i * min_V^b_i  (scalar constant)
                %   h_max = a_i * max_V^b_i  (scalar constant)
                %
                % Constraint 1: p/c >= h_min*u + u_min*a*z - h_min*u_min
                %   <==>  -p/c + h_min*u + u_min*a*z <= h_min*u_min
                %
                % Constraint 2: p/c >= h_max*u + u_max*a*z - h_max*u_max
                %   <==>  -p/c + h_max*u + u_max*a*z <= h_max*u_max
                %
                % Constraint 3: p/c <= h_max*u + u_min*a*z - h_max*u_min
                %   <==>   p/c - h_max*u - u_min*a*z <= -h_max*u_min
                %
                % Constraint 4: p/c <= h_min*u + u_max*a*z - h_min*u_max
                %   <==>   p/c - h_min*u - u_max*a*z <= -h_min*u_max
                % ---------------------------------------------------------
                h_min = s(i).a * s(i).min_V^s(i).b;
                h_max = s(i).a * s(i).max_V^s(i).b;
                u_min = s(i).min_ut;
                u_max = s(i).max_ut;
                a_i   = s(i).a;

                % Constraint 1
                row = row + 1;
                rows  = [rows;  row;             row;            row           ];
                cols  = [cols;  idx_p(i,t)+1;    idx_u(i,t)+1;  idx_z(i,t)+1  ];
                vals  = [vals;  -1/c;            h_min;         u_min*a_i     ];
                rhs   = [rhs;   h_min*u_min       ];
                sense = [sense; '<'];

                % Constraint 2
                row = row + 1;
                rows  = [rows;  row;             row;            row           ];
                cols  = [cols;  idx_p(i,t)+1;    idx_u(i,t)+1;  idx_z(i,t)+1  ];
                vals  = [vals;  -1/c;            h_max;         u_max*a_i     ];
                rhs   = [rhs;   h_max*u_max       ];
                sense = [sense; '<'];

                % Constraint 3
                row = row + 1;
                rows  = [rows;  row;             row;            row           ];
                cols  = [cols;  idx_p(i,t)+1;    idx_u(i,t)+1;  idx_z(i,t)+1  ];
                vals  = [vals;  1/c;            -h_max;        -u_min*a_i     ];
                rhs   = [rhs;  -h_max*u_min       ];
                sense = [sense; '<'];

                % Constraint 4
                row = row + 1;
                rows  = [rows;  row;             row;            row           ];
                cols  = [cols;  idx_p(i,t)+1;    idx_u(i,t)+1;  idx_z(i,t)+1  ];
                vals  = [vals;  1/c;            -h_min;        -u_max*a_i     ];
                rhs   = [rhs;  -h_min*u_max       ];
                sense = [sense; '<'];

                % ---------------------------------------------------------
                % Mass balance: V(i,t) - V(i,t-1) + u(i,t) + sp(i,t) = I(i,t)
                %   [+ u(i-1,t-lag) + sp(i-1,t-lag)  if downstream & t > lag]
                % ---------------------------------------------------------
                if i == 1 || t <= lag
                    row = row + 1;
                    rows  = [rows;  row;            row;             row;            row            ];
                    cols  = [cols;  idx_V(i,t)+1;   idx_V(i,t-1)+1;  idx_u(i,t)+1;  idx_sp(i,t)+1  ];
                    vals  = [vals;  1;              -1;              1;              1              ];
                    rhs   = [rhs;   I(i,t)          ];
                    sense = [sense; '='];
                else
                    row = row + 1;
                    rows  = [rows;  row;            row;             row;            row;            row;                   row                   ];
                    cols  = [cols;  idx_V(i,t)+1;   idx_V(i,t-1)+1;  idx_u(i,t)+1;  idx_sp(i,t)+1; idx_u(i-1,t-lag)+1;   idx_sp(i-1,t-lag)+1   ];
                    vals  = [vals;  1;              -1;              1;              1;              -1;                   -1                    ];
                    rhs   = [rhs;   I(i,t)          ];
                    sense = [sense; '='];
                end

            end  % t == 1 / t >= 2
        end  % unit i
    end  % time t

    nCons = row;
    A = sparse(rows, cols, vals, nCons, nVars);

    % =====================================================================
    %% General Power Constraints: z(i,t) = V(i,t)^b_i  for t = 2..T
    %
    % Struct fields (Gurobi MATLAB API):
    %   .xvar : 0-based index of the argument variable (V)
    %   .yvar : 0-based index of the result variable  (z)
    %   .a    : constant exponent (b_i)
    %
    % IMPORTANT: xvar (V) must have lb >= 0 when b_i is non-integer.
    %            Since min_V > 0 by assumption, this is satisfied.
    % =====================================================================
    genconpow = struct('xvar', {}, 'yvar', {}, 'a', {});
    k = 0;
    for i = 1:n
        for t = 2:T
            k = k + 1;
            genconpow(k).xvar = idx_V(i,t);   % 0-based index of V(i,t)
            genconpow(k).yvar = idx_z(i,t);   % 0-based index of z(i,t)
            genconpow(k).a    = s(i).b;        % exponent in (0,1)
        end
    end

    % =====================================================================
    %% Assemble Gurobi Model
    % =====================================================================
    grb_model.modelsense = 'min';          % Gurobi always minimizes
    grb_model.obj        = obj_coeff;      % -p + 1e-4*sp
    grb_model.A          = A;
    grb_model.rhs        = rhs;
    grb_model.sense      = sense;
    grb_model.lb         = lb;
    grb_model.ub         = ub;
    grb_model.vtype      = repmat('C', nVars, 1);  % all continuous
    if ~isempty(genconpow)
        grb_model.genconpow = genconpow;
    end

    % =====================================================================
    %% Gurobi Parameters
    % FuncNonlinear = 1 : use spatial branch-and-bound (exact, globally
    %                     optimal) instead of static PWL approximation.
    % =====================================================================
    params.OutputFlag    = 1;
    params.Seed          = 1;
    params.Threads       = 1;
    params.FuncNonlinear = 1;

    % =====================================================================
    %% Solve
    % =====================================================================
    result = gurobi(grb_model, params);
    fprintf('Optimization status: %s\n', result.status);

    % =====================================================================
    %% Extract Solution
    % =====================================================================
    if strcmp(result.status, 'OPTIMAL') || strcmp(result.status, 'SUBOPTIMAL')
        x = result.x;

        V_opt  = zeros(n, T);
        p_opt  = zeros(n, T);
        u_opt  = zeros(n, T);
        sp_opt = zeros(n, T);

        for i = 1:n
            for t = 1:T
                V_opt(i,t)  = x(idx_V(i,t)+1);
                p_opt(i,t)  = x(idx_p(i,t)+1);
                u_opt(i,t)  = x(idx_u(i,t)+1);
                sp_opt(i,t) = x(idx_sp(i,t)+1);
            end
        end

        % Reconstruct effective inflow (natural + cascade upstream)
        q = I;
        for i = 2:n
            for t = lag+1:T
                q(i,t) = I(i,t) + u_opt(i-1,t-lag) + sp_opt(i-1,t-lag);
            end
        end

        for i = 1:n
            X = [X, V_opt(i,:)', p_opt(i,:)', u_opt(i,:)', sp_opt(i,:)', q(i,:)'];
        end

        obj = sum(p_opt(:));
    else
        warning('Gurobi did not find an optimal solution. Status: %s', result.status);
        obj = NaN;
    end

end

function [result, obj, X] = realtimeGurobi(t, c, I_t, V_prev, u_prev, V_ref, theta, sys)
% =========================================================================
% INPUTS
%   t        : Current time period index 
%   c        : Power conversion coefficient (scalar)
%   I_t      : Inflow vector for this period (n x 1)
%   V_prev   : Measured storage at start of this period (n x 1) 
%   u_prev   : Turbine release in previous period (n x 1) 
%   V_ref    : Reference storage target from long-term planner (n x 1)
%   theta    : Weight on the volume-tracking penalty term (scalar)
%   sys      : Struct array of per-unit physical parameters:
% =========================================================================

    n = numel(sys);

    % 1: Define Decision Variables
    nVars  = 6 * n;
    idx_V  = @(i)  (i - 1);
    idx_p  = @(i)  n     + (i - 1);
    idx_u  = @(i)  2*n   + (i - 1);
    idx_sp = @(i)  3*n   + (i - 1);
    idx_d  = @(i)  4*n   + (i - 1);
    idx_z  = @(i)  5*n   + (i - 1);   % z(i) = a_i * V(i)^b  (hydraulic head)

    % 2: Define Variable Bounds and Objective Coefficients
    [lb, ub, obj_coeff] = buildVariables(n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_z, V_prev, I_t, theta, sys);

    % 3: Linear Constraints (includes McCormick envelopes)
    [A, rhs, sense] = buildLinearConstraints(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_z, V_prev, u_prev, I_t, V_ref, c, sys);

    % 4: Solve Gurobi Model (pure LP, no nonlinear constraints)
    result = solveModel(n, nVars, obj_coeff, A, rhs, sense, lb, ub, t);

    % 5: Extract Solution
    [obj, X] = extractSolution(result, n, idx_V, idx_p, idx_u, idx_sp, V_prev, I_t, V_ref, c, t, sys);

end


% Helper: Define upper and lower bounds on decision variables & obj. coefs
function [lb, ub, obj_coeff] = buildVariables(n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_z, V_prev, I_t, theta, sys)

    lb        = zeros(nVars, 1);
    ub        =  inf(nVars, 1);
    obj_coeff =  zeros(nVars, 1);

    for i = 1:n

        % B1: Storage V(i)
        lb(idx_V(i)+1)  = sys(i).min_V;
        ub(idx_V(i)+1)  = sys(i).max_V;

        % B2: Power output p(i)
        lb(idx_p(i)+1)  = 0;
        ub(idx_p(i)+1)  = sys(i).F;
        obj_coeff(idx_p(i)+1) = -1;

        % B3: Turbine release u(i)
        lb(idx_u(i)+1)  = sys(i).min_ut;
        ub(idx_u(i)+1)  = sys(i).max_ut;

        % B4: Spill sp(i)
        lb(idx_sp(i)+1) = 0;
        ub(idx_sp(i)+1) = inf;
        obj_coeff(idx_sp(i)+1) = 1e-4;

        % B5: Tracking error d(i)
        lb(idx_d(i)+1)  = 0;
        ub(idx_d(i)+1)  = inf;
        obj_coeff(idx_d(i)+1) = theta / (sys(i).max_V - sys(i).min_V);

        % B6: Hydraulic head z(i) = a_i * V(i)^b
        % Dynamic bounds tightened around V_prev to sharpen McCormick envelope
        V_lo = max(sys(i).min_V, V_prev(i) * 0.8);
        V_hi = min(sys(i).max_V, V_prev(i) * 1.2);
        lb(idx_z(i)+1) = sys(i).a * V_lo^sys(i).b;
        ub(idx_z(i)+1) = sys(i).a * V_hi^sys(i).b;
    end
end


% Helper: Define linear constraints
function [A, rhs, sense] = buildLinearConstraints(t, n, nVars, idx_V, idx_p, idx_u, idx_sp, idx_d, idx_z, V_prev, u_prev, I_t, V_ref, c, sys)

    rows  = [];
    cols  = [];
    vals  = [];
    rhs   = [];
    sense = '';
    row   = 0;

    for i = 1:n

        % (C1) Mass balance: V(i) + u(i) + sp(i) = V_prev(i) + I_t(i)
        row = row + 1;
        rows  = [rows;  row;        row;          row         ];
        cols  = [cols;  idx_V(i)+1; idx_u(i)+1;   idx_sp(i)+1 ];
        vals  = [vals;  1;          1;            1           ];
        rhs   = [rhs;   V_prev(i) + I_t(i)                    ];
        sense = [sense; '='                                    ];

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

        % (C4) Head linearisation: z(i) = a_i * b * V_prev^(b-1) * V(i) + a_i*(1-b)*V_prev^b
        % Tangent of a*V^b at V_prev — exact when V=V_prev, conservative otherwise
        slope = sys(i).a * sys(i).b * V_prev(i)^(sys(i).b - 1);
        const = sys(i).a * (1 - sys(i).b) * V_prev(i)^sys(i).b;
        row = row + 1;
        rows  = [rows;  row;         row        ];
        cols  = [cols;  idx_z(i)+1;  idx_V(i)+1 ];
        vals  = [vals;  1;          -slope      ];
        rhs   = [rhs;   const                   ];
        sense = [sense; '='                     ];

        % (C5-C8) McCormick envelope: p(i) = c * z(i) * u(i)
        % Dynamic bounds tightened around V_prev and ramp window
        V_lo  = max(sys(i).min_V, V_prev(i) * 0.8);
        V_hi  = min(sys(i).max_V, V_prev(i) * 1.2);
        h_min = sys(i).a * V_lo^sys(i).b;
        h_max = sys(i).a * V_hi^sys(i).b;

        if t == 1
            u_lo = sys(i).min_ut;
            u_hi = sys(i).min_ut;   % pinned at t=1
        else
            u_lo = max(sys(i).min_ut, u_prev(i) + sys(i).RR_dn);
            u_hi = min(sys(i).max_ut, u_prev(i) + sys(i).RR_up);
        end

        % MC1: -p/c + h_min*u + u_lo*z <= h_min*u_lo
        row = row + 1;
        rows  = [rows;  row;         row;          row        ];
        cols  = [cols;  idx_p(i)+1;  idx_u(i)+1;   idx_z(i)+1 ];
        vals  = [vals;  -1/c;        h_min;        u_lo       ];
        rhs   = [rhs;   h_min*u_lo                            ];
        sense = [sense; '<'                                    ];

        % MC2: -p/c + h_max*u + u_hi*z <= h_max*u_hi
        row = row + 1;
        rows  = [rows;  row;         row;          row        ];
        cols  = [cols;  idx_p(i)+1;  idx_u(i)+1;   idx_z(i)+1 ];
        vals  = [vals;  -1/c;        h_max;        u_hi       ];
        rhs   = [rhs;   h_max*u_hi                            ];
        sense = [sense; '<'                                    ];

        % MC3: p/c - h_max*u - u_lo*z <= -h_max*u_lo
        row = row + 1;
        rows  = [rows;  row;         row;          row        ];
        cols  = [cols;  idx_p(i)+1;  idx_u(i)+1;   idx_z(i)+1 ];
        vals  = [vals;  1/c;        -h_max;       -u_lo       ];
        rhs   = [rhs;  -h_max*u_lo                            ];
        sense = [sense; '<'                                    ];

        % MC4: p/c - h_min*u - u_hi*z <= -h_min*u_hi
        row = row + 1;
        rows  = [rows;  row;         row;          row        ];
        cols  = [cols;  idx_p(i)+1;  idx_u(i)+1;   idx_z(i)+1 ];
        vals  = [vals;  1/c;        -h_min;       -u_hi       ];
        rhs   = [rhs;  -h_min*u_hi                            ];
        sense = [sense; '<'                                    ];

        % (C9) Tracking error upper:  V(i) - d(i) <= V_ref(i)
        row = row + 1;
        rows  = [rows;  row;         row         ];
        cols  = [cols;  idx_V(i)+1;  idx_d(i)+1  ];
        vals  = [vals;  1;          -1           ];
        rhs   = [rhs;   V_ref(i)                 ];
        sense = [sense; '<'                      ];

        % (C10) Tracking error lower: -V(i) - d(i) <= -V_ref(i)
        row = row + 1;
        rows  = [rows;  row;         row         ];
        cols  = [cols;  idx_V(i)+1;  idx_d(i)+1  ];
        vals  = [vals;  -1;         -1           ];
        rhs   = [rhs;  -V_ref(i)                 ];
        sense = [sense; '<'                      ];

    end

    nCons = row;
    A = sparse(rows, cols, vals, nCons, nVars);
end


% Helper: Solve Gurobi model for the current period (pure LP)
function result = solveModel(n, nVars, obj_coeff, A, rhs, sense, lb, ub, t_period)

    grb_model.modelname  = sprintf('hydroRT_t%d', t_period);
    grb_model.modelsense = 'min';
    grb_model.obj        = obj_coeff;
    grb_model.A          = A;
    grb_model.rhs        = rhs;
    grb_model.sense      = sense;
    grb_model.lb         = lb;
    grb_model.ub         = ub;
    grb_model.vtype      = repmat('C', nVars, 1);

    % Pure LP: no NonConvex or FuncNonlinear needed
    params.OutputFlag = 0;
    params.Seed       = 1;
    params.TimeLimit  = 10;
    params.Threads    = 0;

    result = gurobi(grb_model, params);

    if strcmp(result.status, 'INFEASIBLE')
        iis = gurobi_iis(grb_model);
        iis_rows = find(iis.Arows);
        for r = iis_rows'
            fprintf('Row %d | sense=%s | rhs=%.6f\n', r, grb_model.sense(r), grb_model.rhs(r));
            [~, cols, vals] = find(grb_model.A(r,:));
            for k = 1:length(cols)
                fprintf('  col %d  coeff=%.6f  lb=%.6f  ub=%.6f\n', ...
                    cols(k), vals(k), grb_model.lb(cols(k)), grb_model.ub(cols(k)));
            end
        end
    end

end


% Helper: Extracts the decision variable values
function [obj, X] = extractSolution(result, n, idx_V, idx_p, idx_u, idx_sp, V_prev, I_t, V_ref, c, t_period, sys)

    obj = NaN;
    X   = zeros(n, 5);   % [V, p, u, sp, I] per unit

    if strcmp(result.status, 'OPTIMAL') || strcmp(result.status, 'SUBOPTIMAL')
        x = result.x;

        V_out  = zeros(n, 1);
        p_out  = zeros(n, 1);
        u_out  = zeros(n, 1);
        sp_out = zeros(n, 1);

        for i = 1:n
            V_out(i)  = x(idx_V(i)+1);
            p_out(i)  = x(idx_p(i)+1);
            u_out(i)  = x(idx_u(i)+1);
            sp_out(i) = x(idx_sp(i)+1);
        end

        % Physical power verification
        p_physical = zeros(n, 1);
        for i = 1:n
            h_exact       = sys(i).a * V_out(i)^sys(i).b;
            p_physical(i) = c * h_exact * u_out(i);
        end

        % Normalised volume tracking error
        track_err = zeros(n, 1);
        for i = 1:n
            vol_range    = sys(i).max_V - sys(i).min_V;
            track_err(i) = abs(V_out(i) - V_ref(i)) / vol_range;
        end

        p_error = p_out - p_physical;
        fprintf('[t=%3d]  Power: %7.3f  Spill: %6.3f  PhysErr: %+.2e  TrackErr: %.4f\n', ...
            t_period, sum(p_out), sum(sp_out), max(abs(p_error)), mean(track_err));

        obj = sum(p_out);

        % Store Results: X = [V, p, u, sp, I]
        for i = 1:n
            X(i,:) = [V_out(i), p_out(i), u_out(i), sp_out(i), I_t(i)];
        end

    else
        warning('[t=%d] Gurobi status: %s. Returning NaN obj and zero X.', ...
            t_period, result.status);
    end
end

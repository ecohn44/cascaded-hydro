function [result, obj, X, std_hat, phi_val, alpha_vals] = realtimeGurobi( ...
    t, c, kappa, eps, I_prev, q_error, V_prev, u_prev, V_ref, theta, ...
    lag, up_release, up_ramp, sys, inflow_model, bounds, framework, tracking)

% REALTIMEGUROBINATIVE Native-Gurobi version of realtimeYALMIP.
%
% Variable order in the Gurobi model is
%   z = [V; p; u; sp; d; w].
%
% X retains the original output order
%   X = [V, p, u, sp, q_t].
%

    n = numel(sys);
    bounds    = string(bounds);
    framework = string(framework);

    Rcorr = [1.000000, 0.092792, 0.013465, 0.025220;
             0.092792, 1.000000, 0.159639, 0.090705;
             0.013465, 0.159639, 1.000000, 0.062000;
             0.025220, 0.090705, 0.062000, 1.000000];

    if size(Rcorr,1) ~= n
        error('Rcorr is %d-by-%d, but sys contains %d units.', ...
              size(Rcorr,1), size(Rcorr,2), n);
    end

    %% 1. Forecast inflow and its standard deviation
    q_t = I_prev(:);
    if t > lag
        for i = 2:n
            j = i - 1;  % inflow_model 1=JDA, 2=TDA, 3=BON
            q_t(i) = inflow_model(j).coef0 + inflow_model(j).coef1 * I_prev(i) + inflow_model(j).coef2 * up_release(i-1);
        end
    end

    std_hat = forecast_error(t, kappa, q_error, up_release, up_ramp, ...
                             framework, inflow_model, sys);

    %% 2. Effective storage bounds
    switch bounds
        case {"det", "jcc-ssh"}
            for i = 1:n
                sys(i).V_eff_min = sys(i).min_V;
                sys(i).V_eff_max = sys(i).max_V;
            end

        case "jcc-bon"
            z_score = norminv(1 - eps/(2*n));
            std_V   = [sys.kV]' .* std_hat(:);

            for i = 1:n
                sys(i).V_eff_min = sys(i).min_V + z_score*std_V(i);
                sys(i).V_eff_max = sys(i).max_V;
            end

        otherwise
            error('Unknown bound formulation: %s', bounds);
    end

    %% 3. Build the native Gurobi LP once
    [grb_model, idx] = buildGurobiModel(t, c, V_prev, u_prev, q_t, V_ref, theta, sys);

    params.OutputFlag = 0;
    params.TimeLimit  = 10;
    params.Seed       = 1;
    params.Threads    = 1;

    phi_val    = NaN;
    alpha_vals = zeros(n,1);

    %% 4. Solve deterministic/Bonferroni model or run SSH
    switch bounds
        case {"det", "jcc-bon"}
            [result, obj, X] = solveAndExtractNative( grb_model, idx, q_t, params, false);

        case "jcc-ssh"
            switch framework
                case "det"
                    Sigma_q = zeros(n);

                case {"diu", "ddu"}
                    D_t = diag(std_hat(:));
                    Sigma_q = D_t * Rcorr * D_t;
                    Sigma_q = 0.5 * (Sigma_q + Sigma_q');

                otherwise
                    error('Unknown uncertainty framework: %s', framework);
            end

            if t == 1
                [result, obj, X] = solveAndExtractNative( ...
                    grb_model, idx, q_t, params, false);
                return;
            end

            try
                X_prev_row = zeros(1,5*n);
                for i = 1:n
                    X_prev_row(5*(i-1)+1) = V_prev(i);
                    X_prev_row(5*(i-1)+3) = u_prev(i);
                end

                x_slater = findSlater(X_prev_row, q_t, sys, c);

                [~, x_sol, phi_val, alpha_vals] = applySSHNative( ...
                    grb_model, idx, t, X_prev_row, q_t, Sigma_q, ...
                    x_slater, 1-eps, sys, params);

                result = struct('problem', 0, 'info', 'SSH solved', ...
                                'ssh_failed', false);
                X = compactXToOutput(x_sol, q_t, n);
                obj = sum(X(:,2));

            catch ME
                warning('SSH failed at t=%d: %s. Using deterministic fallback.', ...
                        t, ME.message);

                [fallback, obj, X] = solveAndExtractNative(  grb_model, idx, q_t, params, true);

                if fallback.problem == 0
                    result = struct('problem', 0, 'info', 'Deterministic fallback solved', 'ssh_failed', true);
                else
                    result = struct('problem', 1, 'info', fallback.info, 'ssh_failed', true);
                    X   = [V_prev(:), zeros(n,3), q_t(:)];
                    obj = 0;
                end

                phi_val    = NaN;
                alpha_vals = zeros(n,1);
            end
    end

end


function [m, idx] = buildGurobiModel(t, c, V_prev, u_prev, q_t, V_ref, theta, sys)

    n    = numel(sys);
    nvar = 6*n;

    idx.V  = 1:n;
    idx.p  = n   + (1:n);
    idx.u  = 2*n + (1:n);
    idx.sp = 3*n + (1:n);
    idx.d  = 4*n + (1:n);
    idx.w  = 5*n + (1:n);

    if numel(V_ref) < 2*n
        error('V_ref must contain n lower and n upper reference values.');
    end
    V_lower = V_ref(1:n);
    V_upper = V_ref(n+1:2*n);

    if isscalar(theta)
        theta = repmat(theta, n, 1);
    else
        theta = theta(:);
    end
    if numel(theta) ~= n
        error('theta must be scalar or contain one value per unit.');
    end

    lb = -inf(nvar,1);
    ub =  inf(nvar,1);

    % Eight linear rows per unit: mass balance, four McCormick rows,
    % power equation, and two tracking rows.
    A     = zeros(8*n, nvar);
    rhs   = zeros(8*n, 1);
    sense = repmat(' ', 8*n, 1);
    r = 0;

    for i = 1:n
        VL = sys(i).V_eff_min;
        VU = sys(i).V_eff_max;
        uL = sys(i).min_ut;
        uU = sys(i).max_ut;

        lb(idx.V(i)) = VL;
        ub(idx.V(i)) = VU;
        lb(idx.p(i)) = 0;
        ub(idx.p(i)) = sys(i).F;
        lb(idx.sp(i)) = 0;
        lb(idx.d(i))  = 0;
        lb(idx.w(i))  = 0;
        ub(idx.w(i))  = VU*uU;

        if t == 1
            lb(idx.u(i)) = uL;
            ub(idx.u(i)) = uL;
        elseif u_prev(i) <= 1e-8
            lb(idx.u(i)) = uL;
            ub(idx.u(i)) = uU;
        else
            lb(idx.u(i)) = max(uL, u_prev(i) + sys(i).RR_dn);
            ub(idx.u(i)) = min(uU, u_prev(i) + sys(i).RR_up);
        end

        % Mass balance: V + kV*u + kV*sp = V_prev + kV*q.
        r = r + 1;
        A(r,[idx.V(i), idx.u(i), idx.sp(i)]) = [1, sys(i).kV, sys(i).kV];
        rhs(r)   = V_prev(i) + sys(i).kV*q_t(i);
        sense(r) = '=';

        % McCormick 1: w - VL*u - uL*V >= -VL*uL.
        r = r + 1;
        A(r,[idx.V(i), idx.u(i), idx.w(i)]) = [-uL, -VL, 1];
        rhs(r)   = -VL*uL;
        sense(r) = '>';

        % McCormick 2: w - VU*u - uU*V >= -VU*uU.
        r = r + 1;
        A(r,[idx.V(i), idx.u(i), idx.w(i)]) = [-uU, -VU, 1];
        rhs(r)   = -VU*uU;
        sense(r) = '>';

        % McCormick 3: w - VU*u - uL*V <= -VU*uL.
        r = r + 1;
        A(r,[idx.V(i), idx.u(i), idx.w(i)]) = [-uL, -VU, 1];
        rhs(r)   = -VU*uL;
        sense(r) = '<';

        % McCormick 4: w - VL*u - uU*V <= -VL*uU.
        r = r + 1;
        A(r,[idx.V(i), idx.u(i), idx.w(i)]) = [-uU, -VL, 1];
        rhs(r)   = -VL*uU;
        sense(r) = '<';

        % Power: p - c*b*u - c*a*w = 0.
        r = r + 1;
        A(r,[idx.p(i), idx.u(i), idx.w(i)]) = ...
            [1, -c*sys(i).b, -c*sys(i).a];
        rhs(r)   = 0;
        sense(r) = '=';

        % Tracking: V - d <= V_upper and -V - d <= -V_lower.
        r = r + 1;
        A(r,[idx.V(i), idx.d(i)]) = [1, -1];
        rhs(r)   = V_upper(i);
        sense(r) = '<';

        r = r + 1;
        A(r,[idx.V(i), idx.d(i)]) = [-1, -1];
        rhs(r)   = -V_lower(i);
        sense(r) = '<';
    end

    P_base = sum([sys.F]);
    if P_base <= 0
        error('sum([sys.F]) must be positive.');
    end

    obj = zeros(nvar,1);
    obj(idx.p)  = -1/P_base;
    obj(idx.sp) = 1e-2;
    obj(idx.d)  = theta;

    m.A          = sparse(A(1:r,:));
    m.obj        = obj;
    m.rhs        = rhs(1:r);
    m.sense      = sense(1:r);
    m.lb         = lb;
    m.ub         = ub;
    m.vtype      = repmat('C', nvar, 1);
    m.modelsense = 'min';
    m.modelname  = sprintf('realtime_t%d', t);
end


function [result, obj, X] = solveAndExtractNative(m, idx, q_t, params, ssh_failed)

    grb = gurobi(m, params);
    ok  = isfield(grb,'status') && strcmp(grb.status,'OPTIMAL') ...
          && isfield(grb,'x') && all(isfinite(grb.x));

    if isfield(grb,'status')
        info = grb.status;
    else
        info = 'Gurobi returned no status';
    end

    result = struct('problem', double(~ok), 'info', info, ...
                    'ssh_failed', ssh_failed);
    obj = NaN;
    X   = zeros(numel(q_t),5);

    if ok
        X   = [grb.x(idx.V), grb.x(idx.p), grb.x(idx.u), ...
               grb.x(idx.sp), q_t(:)];
        obj = sum(X(:,2));
    end
end


function X = compactXToOutput(x, q_t, n)
    X = zeros(n,5);
    for i = 1:n
        base = 4*(i-1);
        X(i,:) = [x(base+1), x(base+2), x(base+3), x(base+4), q_t(i)];
    end
end

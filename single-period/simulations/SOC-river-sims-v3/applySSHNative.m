function [model_out, x_sol, phi_val, alpha_vals, dx_hist] = applySSHNative( grb_model, idx, t, X_prev, q_mean, Sigma_q, x_slater, p_target, sys, params)

    % fprintf('  START t=%d\n', t);

    n_units = numel(sys);
    dim_x   = 4*n_units;
    x_slater = x_slater(:);
    q_mean   = q_mean(:)';

    if numel(x_slater) ~= dim_x || any(~isfinite(x_slater))
        error('x_slater must contain 4*n finite values.');
    end

    dx_hist = [];

    % Initial deterministic LP.
    grb = solveOptimal(grb_model, params, 'Initial SSH LP');
    xk  = extractCompactX(grb.x, idx, n_units);

    % Deterministic special case.
    if all(abs(Sigma_q(:)) < 1e-12)
        model_out  = grb_model;
        x_sol      = xk;
        phi_val    = 1;
        alpha_vals = zeros(n_units,1);
        return;
    end

    weights       = zeros(n_units,1);
    ssh_converged = false;
    max_iter      = 20;
    phi_tol       = 5e-4;

    phi_s = computePhiFromX(x_slater, q_mean, Sigma_q, X_prev, sys);
    if phi_s < p_target
        error('Slater point infeasible: phi(x_slater)=%.4f < target=%.4f.', ...
              phi_s, p_target);
    end
    
    % fprintf('Slater point feasible: phi(x_slater)=%.4f >= target=%.4f\n',  phi_s, p_target);

    for iter = 1:max_iter
        phi_k = computePhiFromX(xk, q_mean, Sigma_q, X_prev, sys);

        if phi_k >= p_target - phi_tol
            fprintf('   Converged at iter %d: phi=%.4f\n', iter, phi_k);
            ssh_converged = true;
            break;
        end

       % fprintf('   Iter %d: Unsafe (phi=%.4f). Generating cut...\n', iter, phi_k);

        %% 1. Locate the chance-constraint boundary by bisection
        lam_low = 0;
        lam_high = 1;

        for bis = 1:50
            lam = 0.5*(lam_low + lam_high);
            x_star = (1-lam)*xk + lam*x_slater;
            phi_star = computePhiFromX(x_star, q_mean, Sigma_q, X_prev, sys);

            if abs(phi_star - p_target) < 1e-4
                break;
            elseif phi_star < p_target
                lam_low = lam;
            else
                lam_high = lam;
            end
        end

        %% 2. Finite-difference gradient with respect to r_i = u_i + sp_i
        g_r = zeros(n_units,1);
        for i = 1:n_units
            base = 4*(i-1);
            ju   = base + 3;
            r_i  = x_star(base+3) + x_star(base+4);
            step = 1e-3*max(1,abs(r_i));

            e = zeros(dim_x,1);
            e(ju) = step;

            f_plus  = computePhiFromX(x_star + e, q_mean, Sigma_q, X_prev, sys);
            f_minus = computePhiFromX(x_star - e, q_mean, Sigma_q, X_prev, sys);
            g_r(i)  = (f_plus - f_minus)/(2*step);
        end

        g_scale = max(abs(g_r));
        if g_scale < 1e-10
            error('SSH gradient is numerically zero.');
        end
        g_r = g_r/g_scale;
        weights = abs(g_r);

        %% 3. Append g'*(u+sp) >= g'*(u_star+sp_star)
        cut = sparse(1, size(grb_model.A,2));
        cut(idx.u)  = g_r';
        cut(idx.sp) = g_r';

        cut_rhs = 0;
        for i = 1:n_units
            base = 4*(i-1);
            cut_rhs = cut_rhs + g_r(i)*(x_star(base+3) + x_star(base+4));
        end

        grb_model.A(end+1,:) = cut;
        grb_model.rhs(end+1,1) = cut_rhs;
        grb_model.sense(end+1,1) = '>';

        %% 4. Resolve the LP and use this newly returned solution
        x_old = xk;
        grb = solveOptimal(grb_model, params, sprintf('SSH LP at iteration %d', iter));
        xk = extractCompactX(grb.x, idx, n_units);
        dx_hist(iter,1) = norm(xk - x_old, Inf);
    end

    phi_val = computePhiFromX(xk, q_mean, Sigma_q, X_prev, sys);

    if phi_val < p_target - phi_tol
        warning('SSH did not converge: phi=%.4f. Using feasible Slater point.', phi_val);
        xk = x_slater;
        phi_val = phi_s;
    end

    if ~ssh_converged
        warning('SSH reached max_iter with final phi=%.4f and target=%.4f.', phi_val, p_target);
    end

    if sum(weights) > 0
        alpha_vals = weights/sum(weights);
    else
        alpha_vals = zeros(n_units,1);
    end

    model_out = grb_model;
    x_sol     = xk;
end


function grb = solveOptimal(model, params, label)
    grb = gurobi(model, params);

    if ~isfield(grb,'status')
        error('%s failed: Gurobi returned no status.', label);
    end
    if ~strcmp(grb.status,'OPTIMAL') || ~isfield(grb,'x') ...
            || any(~isfinite(grb.x))
        error('%s failed with Gurobi status %s.', label, grb.status);
    end
end


function x = extractCompactX(z, idx, n)
    x = zeros(4*n,1);
    for i = 1:n
        base = 4*(i-1);
        x(base+1) = z(idx.V(i));
        x(base+2) = z(idx.p(i));
        x(base+3) = z(idx.u(i));
        x(base+4) = z(idx.sp(i));
    end
end


function phi = computePhiFromX(x, q_mean, Sigma_q, X_prev, sys)
    n = numel(sys);
    q_low  = zeros(n,1);
    q_high = zeros(n,1);

    for i = 1:n
        prev_base = 5*(i-1);
        x_base    = 4*(i-1);

        V_prev = X_prev(prev_base+1);
        u_i     = x(x_base+3);
        sp_i    = x(x_base+4);

        q_low(i)  = u_i + sp_i + (sys(i).min_V - V_prev)/sys(i).kV;
        q_high(i) = u_i + sp_i + (sys(i).max_V - V_prev)/sys(i).kV;
    end

    if all(abs(Sigma_q(:)) < 1e-12)
        mu = q_mean(:);
        phi = double(all(mu >= q_low & mu <= q_high));
        return;
    end

    phi = mvncdf(q_low', q_high', q_mean(:)', Sigma_q);
end

function [cons_out, x_sol, phi_val, alpha_vals, dx_hist] = applySSH(cons, vars, t, X_prev, q_mean, Sigma_q, x_slater, p_target, s, Objective, options)

    %fprintf('  START t=%d\n', t);
    n_units = numel(vars.V);
    dim_x   = 4 * n_units;

    % Ensure consistent vector orientation
    x_slater = x_slater(:);
    q_mean   = q_mean(:)';

    % Init convergence tracking
    dx_hist = [];
    
    % Initialize from Solving LP
    sol = optimize(cons, Objective, options); 

    if sol.problem ~= 0
        error('Initial LP solve failed in SSH initialization.');
    end
    
    % Extract starting values from YALMIP objects
    xk = get_current_x(vars, n_units);

    % Return LP solution for deterministic case
    if all(abs(Sigma_q(:)) < 1e-12)
        cons_out  = cons;
        x_sol     = xk;
        phi_val   = 1;                     % deterministic feasibility
        alpha_vals = zeros(n_units,1);     % no risk allocation
        return;
    end

    w = zeros(n_units, 1);
    
    ssh_converged = false;
    iter_count = 0;
    max_iter = 20; % Safety break 
    phi_tol  = 5e-4;

    phi_s = compute_phi_from_x(x_slater, q_mean, Sigma_q, X_prev, s);

    % Check to make sure slater point satifies reliability criteria 
    if phi_s < p_target
        error('Slater point infeasible: φ(x_slater)=%.4f < p_target=%.4f\n', phi_s, p_target);
    %else 
        %fprintf('Slater point feasible: φ(x_slater)=%.4f >= p_target=%.4f\n', phi_s, p_target);
    end 

    while ~ssh_converged && iter_count < max_iter
        iter_count = iter_count + 1;
        
        % Check Reliability of current xk
        phi_k = compute_phi_from_x(xk, q_mean, Sigma_q, X_prev, s);
        
        % Tolerance check 
        if phi_k >= p_target - phi_tol
            fprintf('   Converged at iter %d: φ=%.4f\n', iter_count, phi_k);
            ssh_converged = true;
            break; % EXIT LOOP
        end
        
        %fprintf('   Iter %d: Unsafe (φ=%.4f). Generating cut...\n', iter_count, phi_k);
        
        % Step 1: Bisection Search 
        lam_low = 0; lam_high = 1;
        
        for bis = 1:50
            lam = 0.5 * (lam_low + lam_high);
            
            x_star  = (1 - lam) * xk + lam * x_slater;
            phi_star = compute_phi_from_x(x_star, q_mean, Sigma_q, X_prev, s);
            
            % Stop on the boundary
            if abs(phi_star - p_target) < 1e-4
                break; % Found boundary
            end
            
            if phi_star < p_target
                lam_low = lam; % Too risky, move toward slater
            else
                lam_high = lam; % Too safe, move toward xk
            end
        end
        
        % Step 2a: Calculate Gradient 
        g_u = zeros(n_units,1);
        g_s = zeros(n_units,1);
        g_r = zeros(n_units,1);   % = g_u + g_s

        for i = 1:n_units
            base = 4*(i-1);
            ju = base + 3;
            js = base + 4;
        
            % denominator 
            r_i = x_star(ju) + x_star(js);
            eps_r = 1e-3 * max(1, abs(r_i));
        
            % tolerance
            e = zeros(dim_x,1);
            e(ju) = eps_r;
        
            f1 = compute_phi_from_x( x_star + e, q_mean, Sigma_q, X_prev, s);
        
            f2 = compute_phi_from_x( x_star - e, q_mean, Sigma_q, X_prev, s);
        
            % finite difference
            g_r(i) = (f1 - f2)/(2*eps_r);
        end

        % Improve numerical scaling of the cut
        g_scale = max(abs(g_r));
        
        if g_scale < 1e-10
            error('SSH gradient is numerically zero.');
        end
        
        g_r = g_r/g_scale;

        % Risk attribution weights 
        w = abs(g_r);
                
        % Step 2b: Add the Linear Cut (∇φ'x >= ∇φ'x*)
        lhs = 0;
        rhs = 0;
        
        for i = 1:n_units
            base = 4*(i-1);

            lhs = lhs + g_r(i) * (vars.u(i) + vars.s(i));
            rhs = rhs + g_r(i) * (x_star(base+3) + x_star(base+4));
        end
        
        % Append the new cut
        cons = [cons, lhs >= rhs];
        
        % Step 3: Resolve the LP with the new cut
        x_old = xk;

        optimize(cons, Objective, options);

        if sol.problem ~= 0
            error('SSH LP solve failed at iteration %d.', iter_count);
        end

        % Update xk for the next loop iteration
        xk = get_current_x(vars, n_units);

        % Track convergence metric
        dx_inf = norm(xk - x_old, Inf);
        dx_hist(iter_count) = dx_inf;
    
    end
    
    % Final outputs
    phi_val = compute_phi_from_x(xk, q_mean, Sigma_q, X_prev, s);
    
    if phi_val < p_target - 5e-4
        warning(['SSH did not converge: phi=%.4f. ' ...
                 'Using feasible Slater point.'], phi_val);
    
        xk = x_slater;
        phi_val = phi_s;
    end
    
    cons_out = cons;
    x_sol = xk;

     if ~ssh_converged
        warning('SSH reached max_iter with phi=%.4f and target=%.4f.', phi_val, p_target);
    end
    
    if sum(w) > 0
        alpha_vals = w / sum(w);
    else
        alpha_vals = zeros(n_units, 1);
    end
end


% Helper to extract current numeric values from YALMIP vars
function x_val = get_current_x(vars, n)
    x_val = zeros(4*n, 1);
    for i=1:n
        base = 4*(i-1);
        x_val(base+1) = value(vars.V(i));
        x_val(base+2) = value(vars.p(i));
        x_val(base+3) = value(vars.u(i));
        x_val(base+4) = value(vars.s(i));
    end
end

function phi = compute_phi_from_x(x, q_mean_vec, Sigma_q, X_prev, s)

    n_units = numel(s);

    q_low  = zeros(n_units,1);
    q_high = zeros(n_units,1);

    % Calculate inflow bounds from mass balance equation 
    for i = 1:n_units
        
        % Extract previous volume for unit i
        baseX   = 5*(i-1);
        V_prev  = X_prev(baseX + 1);   % V_i at t-1

        % Extract outflow decision for unit i during initialization 
        base    = 4*(i-1);
        u_i     = x(base + 3);
        s_i     = x(base + 4);

        % Calculate storage-implied inflow band for unit i
        q_low(i)  = u_i + s_i + (s(i).min_V - V_prev)/s(i).kV;
        q_high(i) = u_i + s_i + (s(i).max_V - V_prev)/s(i).kV;
    end

    % disp(Sigma_q)

    % Deterministic special case 
    if all(abs(Sigma_q(:)) < 1e-12)
        % Probability is 1 if mu in [q_low, q_high], else 0
        in_band = all(q_mean_vec >= q_low' & q_mean_vec <= q_high');
        phi = double(in_band);
        return;
    end

    % Multivariate normal probability that q lies in [q_low, q_high]
    phi = mvncdf(q_low', q_high', q_mean_vec, Sigma_q);
end

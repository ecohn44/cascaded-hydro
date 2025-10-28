function [cons_out, x_sol, phi_val, alpha_vals] = applySSH(cons, vars, t, ...
        X_prev, q_mean_vec, Sigma_q, x_slater, p_target, s, Objective, options)

    fprintf('SSH (t=%d)\n', t);

    alpha_vals = [0;0];
    
    %% STEP 1: Solve current deterministic LP (no JCC)
    optimize(cons, -Objective, options);
    
    % Retrieve current solution
    xk = [value(vars.V1); value(vars.p1); value(vars.u1); value(vars.s1);
          value(vars.V2); value(vars.p2); value(vars.u2); value(vars.s2)];
    
    % Evaluate reliability φ(x_k)
    phi_k = compute_phi_from_x(xk, q_mean_vec, Sigma_q, X_prev, s);
    
    % If feasible already (φ ≥ p_target), exit early
    if phi_k >= p_target
        fprintf('Solution xk meets reliability criteria (φ=%.4f ≥ %.4f)\n', phi_k, p_target);
        cons_out = cons; x_sol = xk; phi_val = phi_k;
        return;
    end
    
    %% STEP 2: Interpolation (find boundary point x_star)
    phi_low = phi_k;
    phi_high = compute_phi_from_x(x_slater, q_mean_vec, Sigma_q, X_prev, s);
    
    if phi_high < p_target
        error('Slater point infeasible: φ(x_slater)=%.4f < p_target=%.4f', phi_high, p_target);
    end
    
    lam_low = 0; lam_high = 1;
    for bis = 1:25
        lam = 0.5 * (lam_low + lam_high);
        x_star = (1 - lam) * xk + lam * x_slater;
        phi_star = compute_phi_from_x(x_star, q_mean_vec, Sigma_q, X_prev, s);
        if phi_star < p_target
            lam_low = lam; phi_low = phi_star;
        else
            lam_high = lam; phi_high = phi_star;
        end
        if abs(phi_star - p_target) < 1e-4
            break;
        end
    end
    
    fprintf('   Found boundary point: λ*=%.3f, φ*=%.4f\n', lam, phi_star);
    
    %% STEP 3: Add supporting hyperplane at x_star
    % Compute gradient (finite differences)
    idx = [3 4 7 8]; 
    g = zeros(8,1);
    for i = idx
        eps_i = 1e-6 * max(1, abs(x_star(i))); % adaptive step
        e = zeros(8,1); e(i) = eps_i;
        f1 = compute_phi_from_x(x_star + e, q_mean_vec, Sigma_q, X_prev, s);
        f2 = compute_phi_from_x(x_star - e, q_mean_vec, Sigma_q, X_prev, s);
        g(i) = (f1 - f2) / (2 * eps_i);
    end

    % --- Slack allocation analysis ---
    % Group gradient components by reservoir
    w1 = abs(g(3)) + abs(g(4));   % u1 + s1 sensitivities
    w2 = abs(g(7)) + abs(g(8));   % u2 + s2 sensitivities
    alpha1 = w1 / (w1 + w2);
    alpha2 = w2 / (w1 + w2);
    alpha_vals = [alpha1, alpha2];
    
    fprintf('   Slack allocation: α1=%.2f, α2=%.2f (sum=1)\n', alpha1, alpha2);
    
    if norm(g) < 1e-9
        fprintf('   Gradient ≈ 0 (flat φ region) → skipping cut\n');
    else
        % Build linear form (∇φ(x*)ᵀ (x - x*) ≤ 0)
        lhs = g(1)*vars.V1 + g(2)*vars.p1 + g(3)*vars.u1 + g(4)*vars.s1 + ...
              g(5)*vars.V2 + g(6)*vars.p2 + g(7)*vars.u2 + g(8)*vars.s2;
        rhs = g' * x_star;
    
        % Add supporting hyperplane cut
        cons = [cons, lhs >= rhs];
    
        fprintf('   Added SSH cut: ||∇φ||=%.3e, φ*=%.4f, target=%.4f\n', norm(g), phi_star, p_target);
    end
    
    %% STEP 4: Re-solve LP with the new cut
    optimize(cons, -Objective, options);
    xk = [value(vars.V1); value(vars.p1); value(vars.u1); value(vars.s1);
          value(vars.V2); value(vars.p2); value(vars.u2); value(vars.s2)];
    
    phi_k = compute_phi_from_x(xk, q_mean_vec, Sigma_q, X_prev, s);
    
    fprintf('   After cut: φ(x)=%.4f\n', phi_k);
    
    % Return outputs
    cons_out = cons; 
    x_sol = xk; 
    phi_val = phi_k;
    
end

function phi = compute_phi_from_x(x, q_mean, Sigma_q, X_prev, s)
    u1=x(3); s1=x(4); u2=x(7); s2=x(8);
    q1_low  = s(1).min_V - X_prev(1) + u1 + s1;
    q1_high = s(1).max_V - X_prev(1) + u1 + s1;
    q2_low  = s(2).min_V - X_prev(6) + u2 + s2;
    q2_high = s(2).max_V - X_prev(6) + u2 + s2;
    phi = mvncdf([q1_low; q2_low], [q1_high; q2_high], q_mean, Sigma_q);
end
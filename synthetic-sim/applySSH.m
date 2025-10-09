function [cons_out, x_sol, phi_val] = applySSH(cons, vars, t, ...
        X_prev, q_mean_vec, std_vec, Sigma_q, x_slater, p_target, s, Objective, options)

fprintf('SSH (t=%d)\n', t);

%% Step 1: Initialization 

% Solve the current & deterministic LP (no volume bounds)
optimize(cons, -Objective, options);

% Retrieve solutions 
xk = [value(vars.V1(t)); value(vars.p1(t)); value(vars.u1(t)); value(vars.s1(t));
      value(vars.V2(t)); value(vars.p2(t)); value(vars.u2(t)); value(vars.s2(t))];

% Evaluate joint probability 
phi_k = compute_phi_from_x(xk, q_mean_vec, Sigma_q, X_prev, s);

% Check if solutions meets reliability target 
if phi_k >= p_target
    fprintf('Solution xk meets reliability criteria \n');
    cons_out = cons; x_sol = xk; phi_val = phi_k; return;
end

% Enter while loop
for it = 1:20
    %% Step 2: Interpolation 
    lam = 1;
    while lam > 0
        % Find lambda to interpolate between xs and xk --> return x_star
        x_star = (1-lam)*xk + lam*x_slater;
        if compute_phi_from_x(x_star, q_mean_vec, Sigma_q, X_prev, s) >= p_target, break; end
        lam = lam - 0.1;
    end

    %% Step 3: Add Cut
    % Compute gradient from finite difference method 
    idx = [3 4 7 8]; g = zeros(8,1); eps = 1e-4;
    for i = idx  % (only u1,s1,u2,s2 affect phi)
        e = zeros(8,1); e(i)=eps;
        f1 = compute_phi_from_x(x_star+e,q_mean_vec,std_vec,Sigma_q,X_prev,s);
        f2 = compute_phi_from_x(x_star-e,q_mean_vec,std_vec,Sigma_q,X_prev,s);
        g(i) = (f1-f2)/(2*eps);
    end

    % Add hyperplane to the constraints grad'*(x - x_star) <= 0
    lhs = g(1)*vars.V1(t) + g(2)*vars.p1(t) + g(3)*vars.u1(t) + g(4)*vars.s1(t)+ ...
           g(5)*vars.V2(t) + g(6)*vars.p2(t) + g(7)*vars.u2(t) + g(8)*vars.s2(t);
    cons = [cons, lhs <= g'*x_star];

    %% Step 4: Solve the New LP
    optimize(cons, -Objective, options);
    xk = [value(vars.V1(t)); value(vars.p1(t)); value(vars.u1(t)); value(vars.s1(t));
          value(vars.V2(t)); value(vars.p2(t)); value(vars.u2(t)); value(vars.s2(t))];
    phi_k = compute_phi_from_x(xk, q_mean_vec, std_vec, Sigma_q, X_prev, s);
    if phi_k >= p_target, break; end
end

cons_out = cons; x_sol = xk; phi_val = phi_k;
end


function phi = compute_phi_from_x(x, q_mean, Sigma_q, X_prev, s)
    u1=x(3); s1=x(4); u2=x(7); s2=x(8);
    q1_low  = s(1).min_V - X_prev(1) + u1 + s1;
    q1_high = s(1).max_V - X_prev(1) + u1 + s1;
    q2_low  = s(2).min_V - X_prev(6) + u2 + s2;
    q2_high = s(2).max_V - X_prev(6) + u2 + s2;
    phi = mvncdf([q1_low;q2_low],[q1_high;q2_high],q_mean,Sigma_q);
end
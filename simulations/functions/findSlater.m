function x_slater = findSlater(X_prev, q_mean, s_struct, c)
    % MAKE_SLATER_SUPER_SIMPLE  Very simple Slater point generator
    % Inputs:
    %   X_prev   : previous-state numeric row (X(t-1,:)), length >= 6 (V1 at 1, V2 at 6)
    %   q_mean   : 2x1 vector [q1_mean; q2_mean]
    %   s_struct : structure array s(1), s(2) with fields:
    %              min_V, max_V, min_ut, max_ut, F, a, b (for head/power)
    %   c        : power conversion constant (scalar)
    %
    % Output:
    %   x_slater 8x1 vector [V1; p1; u1; s1; V2; p2; u2; s2]
    
    % Unpack previous volumes and inflow means
    V1_prev = X_prev(1);
    V2_prev = X_prev(6);
    mu1 = q_mean(1); mu2 = q_mean(2);
    
    % Preallocate
    x_slater = zeros(8,1);
    
    % For each unit create midpoint V and midpoint u, then compute spill and p
    for i = 1:2
        if i == 1
            minV = s_struct(1).min_V; maxV = s_struct(1).max_V;
            umin = s_struct(1).min_ut; umax = s_struct(1).max_ut;
            Fcap = s_struct(1).F;
            a = s_struct(1).a; b = s_struct(1).b;
            Vprev = V1_prev; mu = mu1;
        else
            minV = s_struct(2).min_V; maxV = s_struct(2).max_V;
            umin = s_struct(2).min_ut; umax = s_struct(2).max_ut;
            Fcap = s_struct(2).F;
            a = s_struct(2).a; b = s_struct(2).b;
            Vprev = V2_prev; mu = mu2;
        end
    
        % 1) naive midpoint choices
        V_try = 0.5*(minV + maxV);
        u_try = 0.5*(umin + umax);
    
        % 2) compute spill from mass balance: V_prev + mu = V + u + s  -> s = Vprev + mu - V - u
        s_try = Vprev + mu - V_try - u_try;
    
        % 3) if spill negative, reduce u_try (keep it >= umin) to try to make s >= 0
        if s_try < 0
            needed = -s_try;                % how much we must reduce (u+s) to make s zero
            reducible = u_try - umin;       % how much we can reduce u without violating lower bound
            delta_u = min(reducible, needed);
            u_try = u_try - delta_u;
            s_try = Vprev + mu - V_try - u_try; % recompute
        end
    
        % 4) if s still negative (rare), nudge V_try upward (toward maxV) to create room
        if s_try < 0
            % push V_try toward (Vprev + mu - umin) but not past maxV
            V_try = min(maxV, Vprev + mu - umin);
            s_try = Vprev + mu - V_try - u_try;
            % if still negative, set s_try = 0 and clamp u_try to feasible range
            if s_try < 0
                s_try = 0;
                % adjust u so mass-balance holds: u = Vprev + mu - V - s
                u_try = min(max(umin, Vprev + mu - V_try - s_try), umax);
            end
        end
    
        % 5) Ensure V_try strictly inside bounds (add tiny margin)
        eps_margin = 1e-6 * (abs(maxV - minV) + 1);
        V_try = min(max(V_try, minV + eps_margin), maxV - eps_margin);
    
        % 6) Compute head and power; clamp p to feeder capacity
        h_try = a * (V_try^b);                % or call your map_V_to_h if using PWL mapping
        p_try = c * u_try * h_try;
        p_try = min(max(p_try, 0), Fcap);
    
        % 7) store
        if i == 1
            x_slater(1) = V_try; x_slater(2) = p_try; x_slater(3) = u_try; x_slater(4) = s_try;
        else
            x_slater(5) = V_try; x_slater(6) = p_try; x_slater(7) = u_try; x_slater(8) = s_try;
        end
    end

end
function x_slater = findBaseSlater(X_prev, q_mean, s_struct, c)
    % FINDSLATER  Very simple Slater point generator for n units
    %
    % Inputs:
    %   X_prev   : previous-state numeric row (X(t-1,:)),
    %              laid out as [V1 p1 u1 s1 q1  V2 p2 u2 s2 q2  ...]
    %   q_mean   : n×1 vector of inflow means
    %   s_struct : struct array s(i) with fields:
    %              min_V, max_V, min_ut, max_ut, F, a, b
    %   c        : power conversion constant (scalar)
    %
    % Output:
    %   x_slater : 4n×1 vector [V1; p1; u1; s1; V2; p2; u2; s2; ...]
    
    n_units  = numel(s_struct);
    x_slater = zeros(4*n_units, 1);

    for i = 1:n_units
        % Unpack per-unit parameters
        minV = s_struct(i).min_V;
        maxV = s_struct(i).max_V;
        umin = s_struct(i).min_ut;
        umax = s_struct(i).max_ut;
        Fcap = s_struct(i).F;
        a    = s_struct(i).a;
        b    = s_struct(i).b;

        % Previous volume for unit i from X_prev (5 cols per unit)
        baseX = 5*(i-1);
        Vprev = X_prev(baseX + 1);   % V_i at t-1
        mu    = q_mean(i);           % mean inflow for unit i

        % Naive midpoint choices for V and u 
        V_try = 0.5*(minV + maxV);
        u_try = 0.5*(umin + umax);

        % Compute spill from mass balance
        % V_prev + mu = V + u + s  ->  s = Vprev + mu - V - u
        s_try = Vprev + mu - V_try - u_try;

        % If spill negative, reduce u_try down to umin
        if s_try < 0
            needed    = -s_try;            % how much we must reduce (u+s)
            reducible = u_try - umin;      % how much we can reduce u
            delta_u   = min(reducible, needed);
            u_try     = u_try - delta_u;
            s_try     = Vprev + mu - V_try - u_try;  % recompute
        end

        % If s still negative, nudge V_try upward toward maxV
        if s_try < 0
            % push V_try toward (Vprev + mu - umin) but not past maxV
            V_try = min(maxV, Vprev + mu - umin);
            s_try = Vprev + mu - V_try - u_try;

            % if still negative, set s_try = 0 and clamp u_try
            if s_try < 0
                s_try = 0;
                % adjust u so mass balance holds: u = Vprev + mu - V - s
                u_try = min(max(umin, Vprev + mu - V_try - s_try), umax);
            end
        end

        % Ensure V_try strictly inside bounds 
        eps_margin = 1e-6 * (abs(maxV - minV) + 1);
        V_try = min(max(V_try, minV + eps_margin), maxV - eps_margin);

        % Compute head and power; clamp p to feeder capacity 
        h_try = a * (V_try^b);            % or use your rating-curve helper
        p_try = c * u_try * h_try;
        p_try = min(max(p_try, 0), Fcap);

        % Store into x_slater
        baseXsl = 4*(i-1);
        x_slater(baseXsl + 1) = V_try;
        x_slater(baseXsl + 2) = p_try;
        x_slater(baseXsl + 3) = u_try;
        x_slater(baseXsl + 4) = s_try;
    end
end

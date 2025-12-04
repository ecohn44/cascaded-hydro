function x_slater = findSlater(X_prev, q_mean, s_struct, c, X_bon)
% FINDSLATER  Slater point generator for n units, anchored at JCC-BON if available
%
% x_slater = findSlater(X_prev, q_mean, s_struct, c, X_bon)
%
% Inputs:
%   X_prev   : previous-state numeric row (X(t-1,:)),
%              laid out as [V1 p1 u1 s1 q1  V2 p2 u2 s2 q2  ...]
%   q_mean   : n×1 vector of inflow means at time t
%   s_struct : struct array s(i) with fields:
%              min_V, max_V, min_ut, max_ut, F, a, b
%   c        : power conversion constant (scalar)
%   X_bon    : (optional) JCC-BON solution row at time t,
%              same layout as X: [V1 p1 u1 s1 q1 V2 ...]
%
% Output:
%   x_slater : 4n×1 vector [V1; p1; u1; s1; V2; p2; u2; s2; ...]
%
% Logic:
%   - If X_bon is provided, we start from BON's (V, u) for each unit,
%     then adjust slightly to make them strictly interior and satisfy
%     mass balance with (X_prev, q_mean).
%   - If X_bon is not provided, we fall back to a midpoint-based guess.

    n_units  = numel(s_struct);
    x_slater = zeros(4*n_units, 1);

    useBON = (nargin >= 5) && ~isempty(X_bon);

    % Small interior margins for Slater-ness
    epsV = 1e-6;
    epsU = 1e-6;
    epsS = 1e-8;

    for i = 1:n_units
        % Unpack per-unit parameters
        minV = s_struct(i).min_V;
        maxV = s_struct(i).max_V;
        umin = s_struct(i).min_ut;
        umax = s_struct(i).max_ut;
        Fcap = s_struct(i).F;
        a    = s_struct(i).a;
        b    = s_struct(i).b;

        % Indices
        baseX   = 5*(i-1);
        baseXsl = 4*(i-1);

        % Previous volume and mean inflow
        Vprev = X_prev(baseX + 1);   % V_i at t-1
        mu    = q_mean(i);

        % ----- Step 1: Initial (V,u) guess -----
        if useBON
            % Start from BON solution at time t
            V_bon = X_bon(baseX + 1);
            u_bon = X_bon(baseX + 3);

            % Nudge BON values strictly inside bounds
            V_try = min(max(V_bon, minV + epsV*(abs(maxV-minV)+1)), maxV - epsV*(abs(maxV-minV)+1));
            u_try = min(max(u_bon, umin + epsU*(abs(umax-umin)+1)), umax - epsU*(abs(umax-umin)+1));
        else
            % Fallback: midpoint choices for V and u
            V_try = 0.5*(minV + maxV);
            u_try = 0.5*(umin + umax);
        end

        % ----- Step 2: Enforce mass balance & nonnegative spill -----
        % Vprev + mu = V + u + s  =>  s = Vprev + mu - V - u
        s_try = Vprev + mu - V_try - u_try;

        % If spill negative, reduce u_try down toward umin
        if s_try < 0
            needed    = -s_try;           % how much we must reduce (u+s)
            reducible = u_try - umin;     % how much we can reduce u
            delta_u   = min(reducible, needed);

            u_try = u_try - delta_u;
            s_try = Vprev + mu - V_try - u_try;
        end

        % If s still negative, nudge V_try upward toward maxV
        if s_try < 0
            % target V from mass balance with minimal release umin
            V_target = Vprev + mu - umin;

            % keep inside bounds
            V_try = min(max(V_target, minV + epsV*(abs(maxV-minV)+1)), maxV - epsV*(abs(maxV-minV)+1));
            s_try = Vprev + mu - V_try - u_try;

            % if still negative, clamp s to small positive and adjust u
            if s_try < 0
                s_try = epsS;
                u_try = Vprev + mu - V_try - s_try;

                % enforce release bounds strictly
                u_try = min(max(u_try, umin + epsU*(abs(umax-umin)+1)), umax - epsU*(abs(umax-umin)+1));
                % recompute s to keep exact mass balance after clamping u
                s_try = Vprev + mu - V_try - u_try;
                % guarantee non-negative spill
                if s_try < 0
                    s_try = 0;
                end
            end
        end

        % ----- Step 3: Final interior nudges for V,u -----
        DV = abs(maxV - minV) + 1;
        DU = abs(umax - umin) + 1;
        V_try = min(max(V_try, minV + epsV*DV), maxV - epsV*DV);
        u_try = min(max(u_try, umin + epsU*DU), umax - epsU*DU);

        % Re-sync spill with final V,u to satisfy mass balance exactly
        s_try = Vprev + mu - V_try - u_try;
        if s_try < 0
            % If numerical issues: just clip to 0 and slightly shift V
            s_try = 0;
            V_try = Vprev + mu - u_try - s_try;
            V_try = min(max(V_try, minV + epsV*DV), maxV - epsV*DV);
        end

        % Ensure spill is not tiny negative from rounding
        if s_try < 0
            s_try = 0;
        end

        % ----- Step 4: Compute power and cap at F -----
        h_try = a * (V_try^b);
        p_try = c * u_try * h_try;
        p_try = min(max(p_try, 0), Fcap);

        % ----- Store into x_slater -----
        x_slater(baseXsl + 1) = V_try;
        x_slater(baseXsl + 2) = p_try;
        x_slater(baseXsl + 3) = u_try;
        x_slater(baseXsl + 4) = s_try;
    end
end

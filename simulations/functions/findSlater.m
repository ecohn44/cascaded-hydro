function x_slater = findSlater(X_prev, q_t, sys, c)

    n = numel(sys);
    q_t = q_t(:);

    % [V, p, u, spill] for each unit
    x_slater = zeros(4*n, 1);

    for i = 1:n

        % Previous state
        V_prev = X_prev(5*(i-1) + 1);
        u_prev = X_prev(5*(i-1) + 3);

        % Feasible release interval
        u_min = max(sys(i).min_ut, u_prev + sys(i).RR_dn);

        u_max = min(sys(i).max_ut, u_prev + sys(i).RR_up);

        if u_min > u_max
            error('No feasible release interval for unit %d.', i);
        end

        % Curtail to the minimum feasible turbine release
        u = u_min;

        % Maximum storage attainable without spill
        V_no_spill = V_prev + sys(i).kV*(q_t(i) - u);

        % Minimum release cannot prevent lower-bound infeasibility
        if V_no_spill < sys(i).V_eff_min
            error(['No feasible Slater point for unit %d: ' ...
                   'minimum release gives V = %.4f, ' ...
                   'below V_eff_min = %.4f.'], ...
                   i, V_no_spill, sys(i).V_eff_min);
        end

        % Aim for the center of the effective volume range
        V_mid = 0.5*(sys(i).V_eff_min + ...
                     sys(i).V_eff_max);

        % Do not spill if storage is already below the midpoint
        V = min(V_no_spill, V_mid);

        % Spill only the water above the selected interior point
        spill = max(0, ...
            (V_no_spill - V)/sys(i).kV);

        % Power at the Slater point
        p = c*(sys(i).a*V + sys(i).b)*u;

        if p < -1e-6 || p > sys(i).F + 1e-6
            error(['No feasible Slater point for unit %d: ' ...
                   'power %.4f is outside [0, %.4f].'], ...
                   i, p, sys(i).F);
        end

        % Store [V, p, u, spill]
        base = 4*(i-1);
        x_slater(base + (1:4)) = [V; p; u; spill];

    end
end
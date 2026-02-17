function x_slater = findSlater(X_prev, q_mean, sys, c, season)
    
    % findSlater:      Identify a feasible solution to the current LP
    %
    % Inputs:
    %   X_prev   : previous-state numeric row (X(t-1,:)),
    %              laid out as [V1 p1 u1 s1 q1  V2 p2 u2 s2 q2  ...]
    %   q_mean   : n×1 vector of inflow means
    %   s_struct : struct array s(i) with fields:
    %              min_V, max_V, min_ut, max_ut, F, a, b
    %   c        : power conversion constant (scalar)
    %   V_eff    : effective safety volume bounds 
    %
    % Output:
    %   x_slater : 4n×1 vector [V1; p1; u1; s1; V2; p2; u2; s2; ...]
    
    n_units  = numel(sys);
    x_slater = zeros(4*n_units, 1);
    
    for i = 1:n_units
        
        % Unpack per-unit parameters
        minV = sys(i).min_V; 
        maxV = sys(i).max_V;
        umin = sys(i).min_ut;
        umax = sys(i).max_ut;
        Fcap = sys(i).F;
        a    = sys(i).a;
        b    = sys(i).b;
        RRdn = sys(i).RR_dn;

        % Previous volume for unit i (V_i at t-1)
        baseX = 5*(i-1);
        Vprev = X_prev(baseX + 1);   
       
        % mean inflow for unit i
        mu    = q_mean(i); 

        % Previous water release for unit i (u_i at t-1)
        Uprev = X_prev(baseX + 3);

        if season == "dry"
            % Calculate current volume before release
            V_pre = Vprev + mu;
         
            % Conservative release is full curtailment 
            u_try = max(umin, Uprev + RRdn);
    
            % Calculate spill if volume is above upper bounds
            s_try = 0;
            V_after = V_pre - u_try;
    
            if V_after > maxV
                s_try = V_after - maxV;
                V_after = maxV;
            end 
       
            % Compute head and power; clamp p to feeder capacity 
            V_try = V_after;
            h_try = a * (V_after^b);          
            p_try = c * u_try * h_try;
            p_try = min(max(p_try, 0), Fcap);
        
        elseif season == "wet"
            % Pick target volume total outflow r = u+s ----
            Vmid = 0.5*(minV + maxV);

            % Calculate total outflow r to mean volume target
            r_center = Vprev + mu - Vmid;

            % Feasible r range to keep next volume in [minV,maxV]
            r_lo = Vprev + mu - maxV;
            r_hi = Vprev + mu - minV;

            % Clamp release into feasible interval
            r_try = min(max(r_center, r_lo), r_hi);

            % Partition release into spill and generation 
            u_try = min(max(umin, Uprev + RRdn), r_try);     % propose conservartive turbine release
            s_try = max(0, r_try - u_try);                   % spill absorbs the remainder

            % Calculate resulting volume state 
            V_try = Vprev + mu - (u_try + s_try);

            % Calculate power generation 
            h_try = a * (V_try^b);          
            p_try = c * u_try * h_try;

        end


        % Store into x_slater
        baseXsl = 4*(i-1);
        x_slater(baseXsl + 1) = V_try;
        x_slater(baseXsl + 2) = p_try;
        x_slater(baseXsl + 3) = u_try;
        x_slater(baseXsl + 4) = s_try;

    end
end

%% Policy Test Framework 
% =======================================================================
% Author: Eliza Cohn
% Description: Simulate the reservoir volume trajecotry under historical
%               inflow and optimal release policies learned during training
%               

function [V_sim, u_sim, p_sim, IVI] = runPolicyTestSims(sysparams, bounds, X, policyLabel)

    % Dimensions
    T     = size(X,1);              % number of time steps
    n     = numel(sysparams);       % number of units
    tt    = (1:T)';

    % Extract deterministic (solved) policy for all units
    % X columns per unit i: [V_i, p_i, u_i, s_i, q_mean_i]
    V_opt = zeros(T, n);
    u_opt = zeros(T, n);
    s_opt = zeros(T, n);
    q_hist = zeros(T, n);

    for i = 1:n
        base = 5*(i-1);
        V_opt(:, i)   = X(:, base+1);
        u_opt(:, i)   = X(:, base+3);
        s_opt(:, i)   = X(:, base+4);

        % Extract historical streamflow data 
        q_hist(:, i)  = X(:, base+5);
    end

    % Extract initial volumes and bounds
    Vmin = zeros(1,n);
    Vmax = zeros(1,n);
    V0   = zeros(1,n);
    a    = zeros(1,n);
    b    = zeros(1,n);

    for i = 1:n
        Vmin(i) = sysparams(i).min_V;
        Vmax(i) = sysparams(i).max_V;
        V0(i)   = sysparams(i).V0;
        a(i)    = sysparams(i).a;
        b(i)    = sysparams(i).b;
    end

    % Allocate arrays for historical policy testing
    V_sim = zeros(T, n);
    u_sim = zeros(T, n);
    p_sim = zeros(T, n);
    s_sim = zeros(T, n);

    % Violation indicators: T x n
    violation_count = false(T, n);
    viol_mag = zeros(T, n);

    %% Begin Simulation 
    % Previous volumes for each unit
    V_prev = V0;

    for t = 1:T
        for i = 1:n
            
            % Record historical inflow and sample forecast uncertainty
            inflow_i = q_hist(t,i);

            % Optimal release policy at this time
            u_i = u_opt(t,i);

            % Optimal spill policy 
            % s_i = s_opt(t,i);

            % Calculate spill 
            s_i = max(0, V_prev(i) - Vmax(i) + inflow_i - u_i);

            % Next volume state
            V_i = V_prev(i) + inflow_i - u_i - s_i;

            % Record violations
            below = max(Vmin(i) - V_i, 0);
            above = max(V_i - Vmax(i), 0);
            
            viol_mag(t, i) = below + above;
            violation_count(t, i) = (viol_mag(t, i) > 0);

            if V_i <= 0
                p_i = 0; 
            else
                % Compute head and generation
                h_i = a(i) * (V_i^b(i));
                p_i = u_i * h_i;    % (TEMP) multiply by c when rescaled
            end

            % Store trajectories
            V_sim(t, i) = V_i;
            u_sim(t, i) = u_i;
            p_sim(t, i) = p_i;
            s_sim(t, i) = s_i;

            % Advance state
            V_prev(i) = V_i;
        end
    end

    %% Flood Metrics
    headroom = zeros(T,n);
    for i = 1:n
        headroom(:,i) = Vmax(i) - V_sim(:,i);
    end
    
    AFV_tot  = sum(headroom, 1);   % 1 x n  (sum over time)
    AFV_mean = AFV_tot(:);         % n x 1
    AFV_sys  = sum(AFV_tot);       % scalar
    
    fprintf('\n  AFV (time-integrated headroom) Report\n');
    for i = 1:n
        fprintf('Reservoir %d: AFV = %.4f\n', i, AFV_mean(i));
    end
    fprintf('System Total: AFV = %.4f\n', AFV_sys);


    %% General Simulation Metrics

    % Terminal volume (end-of-horizon storage)
    V_T = V_sim(end, :);   % 1 x n
    
    fprintf('\n  Terminal Storage Report\n');
    for i = 1:n
        fprintf('Reservoir %d: V_T = %.4f\n', i, V_T(i));
    end
    fprintf('System Mean: mean(V_T) = %.4f\n', mean(V_T));
    

    % Metric (IVI): Integrated (time-avg) violation magnitude
    IVI = mean(viol_mag, 1); 

    fprintf('\n  Monte Carlo Violation Report\n');
    
    for i = 1:n
        
        fprintf('Reservoir %d:\n', i);
        fprintf('   Integrated Violation Index (IVI): %.4f\n', IVI(i));
    end
    
    switch string(bounds)
        case "det",     bLabel = "Deterministic";
        case "icc",     bLabel = "Individual CC";
        case "jcc-bon", bLabel = "Bonferroni JCC";
        otherwise,      bLabel = char(bounds);
    end

    % Plot policy test volume trajectories for each reservoir
    % plotPolicyTestVolumes(tt, V_opt, V_sim, Vmin, Vmax, bLabel, policyLabel);
    % plotSpill(tt, s_opt, s_sim, bLabel, policyLabel)

end


function plotPolicyTestVolumes(tt, V_opt, V_sim, Vmin, Vmax, bLabel, policyLabel)

    [~, n] = size(V_sim);

    figure('Name', sprintf('Policy Replay Storage Trajectories — %s', bLabel), ...
       'NumberTitle','off', ...
       'Units','normalized', ...
       'Position',[0.1 0.3 0.8 0.4]);   % [left bottom width height]
    tiledlayout(1, n, 'TileSpacing','compact', 'Padding','compact');


    for i = 1:n
        nexttile; hold on; grid on;
    
        Vp = V_opt(:,i);
        Vs = V_sim(:,i);
    
        % Plot trajectories
        plot(tt, Vp, 'r-', 'LineWidth', 2.0, 'DisplayName','Planned (training)');
        plot(tt, Vs, 'b-', 'LineWidth', 2.0, 'DisplayName','Replay (historical)');
    
        % IVI shading (robust to multiple violation segments)
        idxU = Vs > Vmax(i);
        
        d = diff([false; idxU; false]);
        starts = find(d == 1);
        ends   = find(d == -1) - 1;
        
        for k = 1:numel(starts)
            ii = starts(k):ends(k);
        
            xU = tt(ii);
            y1 = Vmax(i) * ones(size(xU));
            y2 = Vs(ii);
        
            hIVI = fill([xU; flipud(xU)], [y1; flipud(y2)], ...
                [0.6 0.6 0.6], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        
            if i == 1 && k == 1
                set(hIVI, 'DisplayName', 'IVI (V > V_{max})');
            else
                set(hIVI, 'HandleVisibility', 'off');
            end
        end

        % Bounds
        yline(Vmin(i), '--k', 'LineWidth', 1.5, 'HandleVisibility','off');
        yline(Vmax(i), '--k', 'LineWidth', 1.5, 'HandleVisibility','off');

    
        xlabel('Hour');
        ylabel(sprintf('V_%d (p.u.)', i));
        title(sprintf('Reservoir %d', i));
    
        % Limits (pad a bit)
        y_low  = min([Vp; Vs; Vmin(i)]) - 0.03;
        y_high = max([Vp; Vs; Vmax(i)]) + 0.03;
        ylim([y_low y_high]);
        xlim([1, length(tt)])
    
        if i == 1
            legend('Location','best');
        end
    end

    sgtitle(sprintf('Policy Replay Storage Trajectories — %s — %s', policyLabel, bLabel));
end



function plotSpill(tt, s_policy, s_test, bLabel, policyLabel)

    [~, n] = size(s_policy);

    figure('Name', sprintf('Spill Trajectories — %s', bLabel), ...
       'NumberTitle','off', ...
       'Units','normalized', ...
       'Position',[0.1 0.3 0.8 0.4]);   % [left bottom width height]
    tiledlayout(1, n, 'TileSpacing','compact', 'Padding','compact');


    for i = 1:n
        nexttile; hold on; grid on;
    
        Sp = s_policy(:,i);
        St = s_test(:,i);
    
        % Plot trajectories
        plot(tt, Sp, 'r-', 'LineWidth', 2.0, 'DisplayName','Planned (training)');
        plot(tt, St, 'b-', 'LineWidth', 2.0, 'DisplayName','Replay (historical)');
    
        xlabel('Hour');
        ylabel(sprintf('V_%d (p.u.)', i));
        title(sprintf('Reservoir %d', i));
    
        xlim([1, length(tt)])
    
        if i == 1
            legend('Location','best');
        end
    end

    sgtitle(sprintf('Spill Trajectories — %s — %s', policyLabel, bLabel));
end


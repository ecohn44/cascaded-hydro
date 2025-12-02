function [V_sim, u_sim, p_sim, MFV, RLR] = runMonteCarloSims(sysparams, bounds, std_hat, X, savePath, printplot, policyLabel)

    % Set random seed
    rng(0, 'twister');

    % Dimensions
    T     = size(X,1);              % number of time steps
    n     = numel(sysparams);       % number of units
    nSim  = 5000;                   % number of Monte Carlo runs
    tt    = (1:T)';

    % Extract deterministic (solved) policy for all units
    % X columns per unit i: [V_i, p_i, u_i, s_i, q_mean_i]
    V_opt = zeros(T, n);
    u_opt = zeros(T, n);
    s_opt = zeros(T, n);
    q_mean = zeros(T, n);

    for i = 1:n
        base = 5*(i-1);
        V_opt(:, i)   = X(:, base+1);
        u_opt(:, i)   = X(:, base+3);
        s_opt(:, i)   = X(:, base+4);
        q_mean(:, i)  = X(:, base+5);
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

    % Allocate arrays for MC trajectories
    V_sim = zeros(T, n, nSim);
    u_sim = zeros(T, n, nSim);
    p_sim = zeros(T, n, nSim);

    % Violation indicators: T x nSim x n
    violation_count = false(T, nSim, n);

    % MAIN MONTE CARLO LOOP
    for sIdx = 1:nSim
        
        % Previous volumes for each unit
        V_prev = V0;

        for t = 1:T
            for i = 1:n
                % Sample inflow for unit i
                inflow_i = q_mean(t,i) + std_hat(t,i)*randn;

                % Fixed (solved) policy at this time
                u_i = u_opt(t,i);
                s_i = s_opt(t,i);

                % Next volume 
                V_i = V_prev(i) + inflow_i - u_i - s_i;

                % Record violation based on unclamped state
                violation_count(t, sIdx, i) = (V_i < Vmin(i)) || (V_i > Vmax(i));

                if V_i <= 0
                    p_i = 0; 
                else
                    % Compute head and generation
                    h_i = a(i) * (V_i^b(i));
                    p_i = u_i * h_i;   % (TEMP) multiply by c when rescaled
                end

                % Store trajectories
                V_sim(t, i, sIdx) = V_i;
                u_sim(t, i, sIdx) = u_i;
                p_sim(t, i, sIdx) = p_i;

                % Advance state
                V_prev(i) = V_i;
            end
        end
    end

    %% Simulation Metrics
    % Metric #1 (MFV): Per-time-step violation fraction
    viol_frac = squeeze(mean(violation_count, 2));
    MFV = mean(viol_frac, 1);                
    
    % Metric #2 (RLR): Run-level violation probability
    run_violation = squeeze(any(violation_count, 1));  
    RLR = mean(run_violation, 1);     
    
    fprintf('\n  Monte Carlo Violation Report\n');
    
    for i = 1:n
        
        fprintf('Reservoir %d:\n', i);
        
        fprintf('   Mean Frequency of Violations (MFV): %.2f%%\n', ...
                100 * MFV(i));

        fprintf('   Run-Level Risk (RLR): %.2f%%\n', ...
                100 * RLR(i));
    end


    switch string(bounds)
        case "det",     bLabel = "Deterministic";
        case "icc",     bLabel = "Individual CC";
        case "jcc-bon", bLabel = "Bonferroni JCC";
        otherwise,      bLabel = char(bounds);
    end

    % Plot MC volume trajectories for each reservoir
    if printplot
        plotReservoirMCGrid(tt, V_opt, V_sim, Vmin, Vmax, ...
                            bLabel, policyLabel, savePath, bounds, printplot);
    end
end


function plotReservoirMCGrid(tt, V_opt, V_sim, Vmin, Vmax, ...
                             bLabel, policyLabel, savePath, bounds, printplot)

    [~, n, ~] = size(V_sim);

    figure('Name', sprintf('MC Volumes — %s', bLabel), ...
           'NumberTitle', 'off');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:n
        nexttile; hold on; grid on;

        % T x nSim slice for this unit
        V_sim_i = squeeze(V_sim(:, i, :));

        % Mean trajectory
        V_mean = mean(V_sim_i, 2);

        % IQR band
        pr1     = prctile(V_sim_i, [25 75], 2);
        V_lower = pr1(:,1);
        V_upper = pr1(:,2);

        % IQR band fill
        fill([tt; flipud(tt)], [V_upper; flipud(V_lower)], ...
             [0.3 0.5 1.0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        % Optimal and MC mean
        plot(tt, V_opt(:,i), 'r-',  'LineWidth', 2,   'DisplayName','Optimal');
        plot(tt, V_mean,     'b--', 'LineWidth', 1.5, 'DisplayName','MC mean');

        % Bounds
        yline(Vmin(i), '--r');
        yline(Vmax(i), '--r');

        xlabel('Hour');
        ylabel(sprintf('V_%d (p.u.)', i));
        title(sprintf('Reservoir %d', i));

        if i == 1
            legend('IQR band','Optimal','MC mean','Location','best');
        end
    end

    sgtitle(sprintf('Monte Carlo Storage Trajectories — %s — %s', ...
            policyLabel, bLabel));

    if printplot
        fname = sprintf('mc_grid_%s_%s.png', bLabel, char(bounds));
        saveas(gcf, fullfile(savePath, fname));
    end
end

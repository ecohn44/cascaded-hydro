function [V1_sim, V2_sim, u1_sim, u2_sim, p1_mean, p2_mean, deficits1, deficits2] = runMonteCarloSims(sysparams, bounds, std_hat, X, savePath, printplot)
% runMonteCarloSims - Monte Carlo bound violation checker (FIXED POLICY, CLAMPED, WITH SIM HEAD)
%
%   Changes from original version:
%     (1) Added clamping to release/spill so volumes stay within bounds 
%     (2) Added generation recalculation fmor adjusted water release
%     (3) Added average generation reporting across runs
%
%   Inputs:
%       sysparams : struct array (size 2) with fields .min_V, .max_V, .V0, .a, .b
%       bounds    : label for the bound type (e.g., 'det', 'icc', 'jcc-bon')
%       std_hat   : Tx2 matrix of inflow std deviations (same units as q)
%       X         : Tx10 matrix of deterministic results
%                   Columns:
%                     1=V1, 2=p1, 3=u1, 4=s1, 5=q1,
%                     6=V2, 7=p2, 8=u2, 9=s2, 10=q2
%       savePath  : folder for plots
%       printplot : true/false to save figures
%
%   Outputs:
%       V1, V2    : TxN simulated (CLAMPED) volume trajectories

    %  Simulation parameters 
    nSim = 5000;               % number of Monte Carlo runs
    T = size(X,1);             % number of time steps
    tt = (1:T)';               % time vector

    %  Extract system parameters 
    Vmin1 = sysparams(1).min_V; Vmax1 = sysparams(1).max_V; V0_1 = sysparams(1).V0;
    Vmin2 = sysparams(2).min_V; Vmax2 = sysparams(2).max_V; V0_2 = sysparams(2).V0;
    Umin1 = 0; % sysparams(1).min_ut;
    Umin2 = 0; % sysparams(2).min_ut; 
    Umax1 = sysparams(1).max_ut; Umax2 = sysparams(2).max_ut; % will use for wet season
    s_min1 = 0; s_min2 = 0;

    %  Extract deterministic (solved) policy 
    V1_opt = X(:,1);  p1_opt = X(:,2);  u1_opt = X(:,3);  s1_opt = X(:,4);  q1_mean = X(:,5);
    V2_opt = X(:,6);  p2_opt = X(:,7);  u2_opt = X(:,8);  s2_opt = X(:,9);  q2_mean = X(:,10);

    %  Extract inflow uncertainty (std) 
    std1 = std_hat(:,1);
    std2 = std_hat(:,2);

    %  Define head functions from sysparams
    predict_h1 = @(V) sysparams(1).a .* (V.^sysparams(1).b);
    predict_h2 = @(V) sysparams(2).a .* (V.^sysparams(2).b);

    %  Allocate arrays for MC policies (possible curtailment)
    V1_sim = zeros(T, nSim);
    V2_sim = zeros(T, nSim);
    p1_sim = zeros(T, nSim);
    p2_sim = zeros(T, nSim);
    u1_sim = zeros(T, nSim);
    u2_sim = zeros(T, nSim);

    % Original violation counts from unaltered policies 
    violation_count1 = false(T, nSim);
    violation_count2 = false(T, nSim);

    % Curtailed violation counts from saturated policies
    violation_curtail_count1 = false(T, nSim);
    deficits1 = zeros(T, nSim);
    violation_curtail_count2 = false(T, nSim);
    deficits2 = zeros(T, nSim);

    % ===================================================================
    %  MAIN MONTE CARLO LOOP
    % ===================================================================
    for n = 1:nSim
        V1_prev = V0_1; V2_prev = V0_2;

        for t = 1:T
            % Sample inflows
            inflow1 = q1_mean(t) + std1(t)*randn; % To do: DDU needs to update based on updated policy 
            inflow2 = q2_mean(t) + std2(t)*randn;

            % Fixed (solved) policy
            u1 = u1_opt(t); s1 = s1_opt(t);
            u2 = u2_opt(t); s2 = s2_opt(t);

            %  Compute next state update for auditing 
            V1_unc = V1_prev + inflow1 - u1 - s1;
            V2_unc = V2_prev + inflow2 - u2 - s2;

            % Record violations based on UNCLAMPED state
            violation_count1(t,n) = (V1_unc < Vmin1) || (V1_unc > Vmax1);
            violation_count2(t,n) = (V2_unc < Vmin2) || (V2_unc > Vmax2);

            % ===================================================================
            %  State regulation via control saturation 
            % ===================================================================
            [V1, u1, deficits1(t,n), violation_curtail_count1(t,n)] = curtail_deficit(V1_unc, u1_opt(t), Vmin1, Umin1);
            [V2, u2, deficits2(t,n), violation_curtail_count2(t,n)] = curtail_deficit(V2_unc, u2_opt(t), Vmin2, Umin2);


            % ===================================================================
            %  Recalculate generation using simulated head
            % ===================================================================
            h1 = predict_h1(V1);
            h2 = predict_h2(V2);
            p1_sim(t,n) = u1 * h1; % (TEMP) Add back in c when full scale 
            p2_sim(t,n) = u2 * h2;

            % Store curtailed release and clamped volumes
            V1_sim(t,n) = V1; V2_sim(t,n) = V2;
            u1_sim(t,n) = u1; u2_sim(t,n) = u2;

            % Advance to next state
            V1_prev = V1; V2_prev = V2;
        end
    end

    % ===================================================================
    %  Report violation rates and generation stats
    % ===================================================================

    %  Violation audit (based on UNCLAMPED) 
    avg_per_timestep1 = mean(violation_count1, 2);
    avg_per_timestep2 = mean(violation_count2, 2);
    overall_avg1 = mean(avg_per_timestep1);
    overall_avg2 = mean(avg_per_timestep2);

    fprintf('\n  Monte Carlo Violation Report (UNCLAMPED audit) \n');
    fprintf('Reservoir 1: %.2f%%\n', 100*overall_avg1);
    fprintf('Reservoir 2: %.2f%%\n', 100*overall_avg2);

    %  Violation audit (based on CURTAILED)
    c_avg_per_timestep1 = mean(violation_curtail_count1, 2);
    c_avg_per_timestep2 = mean(violation_curtail_count2, 2);
    c_overall_avg1 = mean(c_avg_per_timestep1);
    c_overall_avg2 = mean(c_avg_per_timestep2);

    fprintf('\n  Monte Carlo Violation Report (CURTAILED audit) \n');
    fprintf('Reservoir 1: %.2f%%\n', 100*c_overall_avg1);
    fprintf('Reservoir 2: %.2f%%\n', 100*c_overall_avg2);

    % Compute mean total generation across all runs 
    mean_p1_MC = mean(sum(p1_sim,1));
    mean_p2_MC = mean(sum(p2_sim,1));
    total_p_MC = mean_p1_MC + mean_p2_MC;

    % Compute average generation profile for Unit 02
    p2_mean = mean(p2_sim,2);
    p1_mean = mean(p1_sim,2);

    % Compare with deterministic baseline 
    sum_p1_opt = sum(p1_opt);  sum_p2_opt = sum(p2_opt);
    total_p_opt = sum_p1_opt + sum_p2_opt;

    % Calculate change in power production 
    drop_p1 = 100*(sum_p1_opt - mean_p1_MC)/max(sum_p1_opt,eps);
    drop_p2 = 100*(sum_p2_opt - mean_p2_MC)/max(sum_p2_opt,eps);
    total_change = 100*(total_p_MC - total_p_opt)/total_p_opt;

    fprintf('\n  Average Generation from Simulated Head (CLAMPED) \n');
    fprintf('Res1: mean p = %.5f (vs %.5f opt) | Δ=%.2f%%\n', mean_p1_MC, sum_p1_opt, drop_p1);
    fprintf('Res2: mean p = %.5f (vs %.5f opt) | Δ=%.2f%%\n', mean_p2_MC, sum_p2_opt, drop_p2);
    fprintf('Total Generation (opt): %.3f\n', total_p_opt)
    fprintf('Total Generation (curtail): %.3f\n',total_p_MC)
    fprintf('Percent change in Generation: %.2f\n', total_change)

    stats1 = analyzeFlowMinViolations(u1_sim, sysparams(1).min_ut, "Unit 1");
    stats2 = analyzeFlowMinViolations(u2_sim, sysparams(2).min_ut, "Unit 2");

    % ===================================================================
    %  Visualization of Monte Carlo Runs
    % ===================================================================
    V1_mean = mean(V1_sim,2);  V1_std = std(V1_sim,0,2);
    V2_mean = mean(V2_sim,2);  V2_std = std(V2_sim,0,2);

    % Calculate 1-sigma enevlope accross runs for Volume trajectory 
    V1_upper = V1_mean + V1_std; V1_lower = V1_mean - V1_std;
    V2_upper = V2_mean + V2_std; V2_lower = V2_mean - V2_std;

    switch string(bounds)
        case "det",     bLabel = "Deterministic";
        case "icc",     bLabel = "Individual CC";
        case "jcc-bon", bLabel = "Bonferroni JCC";
        otherwise,      bLabel = char(bounds);
    end

    %  Plot Reservoir 1 
    figure; hold on; grid on;
    fill([tt; flipud(tt)], [V1_upper; flipud(V1_lower)], [0.3 0.5 1.0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(tt, V1_opt, 'r-', 'LineWidth', 2, 'DisplayName','Optimal V1');
    plot(tt, V1_mean, 'b--', 'LineWidth', 1.5, 'DisplayName','MC mean');
    yline(Vmin1, '--r', 'Min bound'); yline(Vmax1, '--r', 'Max bound');
    xlabel('Hour'); ylabel('Storage V_1 (p.u.)');
    title(sprintf('Reservoir 1 — %s', bLabel));
    legend('1\sigma band','Optimal','MC mean','Location','best');
    if printplot, saveas(gcf, fullfile(savePath, ['mc_reservoir1_' char(bounds) '_fixed_clamped.png'])); end

    %  Plot Reservoir 2 
    figure; hold on; grid on;
    fill([tt; flipud(tt)], [V2_upper; flipud(V2_lower)], [0.3 0.5 1.0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(tt, V2_opt, 'r-', 'LineWidth', 2, 'DisplayName','Optimal V2');
    plot(tt, V2_mean, 'b--', 'LineWidth', 1.5, 'DisplayName','MC mean');
    yline(Vmin2, '--r', 'Min bound'); yline(Vmax2, '--r', 'Max bound');
    xlabel('Hour'); ylabel('Storage V_2 (p.u.)');
    title(sprintf('Reservoir 2 — %s', bLabel));
    legend('1\sigma band','Optimal','MC mean','Location','best');
    if printplot, saveas(gcf, fullfile(savePath, ['mc_reservoir2_' char(bounds) '_fixed_clamped.png'])); end

    % ===================================================================
    %  Visualization of  Curtailment & Umin Violations
    % ===================================================================

    plotCurtailment(deficits1, u1_opt, sysparams(1).min_ut, 'Unit 1 Curtailment');
    plotCurtailment(deficits2, u2_opt, sysparams(2).min_ut, 'Unit 2 Curtailment');

end


function [V, u, take_turb, viol] = curtail_deficit(V_unc, u_opt, Vmin, u_min)

    % Record current volume and release
    V = V_unc; u = u_opt; take_turb = 0; viol = 0; 

    % Check if original release policy violated volume bounds
    if V_unc < Vmin
        
        % Calculate deviation from lower bounds
        deficit = Vmin - V_unc;
    
        % Check if release needs to be curtailed
        if deficit > 0
            % Calculate maximum possible curtailment 
            curtail_cap = max(0, u_opt - u_min);
            
            % Compare deficit to max curtailment
            take_turb = min(deficit, curtail_cap);
            
            % Saturate control
            if take_turb > 0
                % Reduce u1
                u = u_opt - take_turb; 
                
                % Update volume to represent curtailed release
                V = V_unc + take_turb;
                
                % Update deficit 
                deficit_remain = deficit - take_turb;
            end
        end
    
        % If deficit still exists, then violation is unavoidable
        if deficit_remain > 0
            viol = 1; 
        end
    end
end


function plotCurtailment(deficits, u_opt, u_min, titleStr)
% Bars = mean(deficits); upper-only whisker = (p99 - mean)
% Red shade = violation zone y >= max(0, u_opt - u_min)

    if nargin < 4, titleStr = 'Curtailment (mean + 99th percentile upper whisker)'; end

    T   = size(deficits,1);
    tt  = (1:T)';

    % stats
    mu = mean(deficits, 2);

    % 99th percentile along simulations (dim 2); fallback if prctile missing
    if exist('prctile','file')
        p99 = prctile(deficits, 99, 2);
    else
        % simple fallback
        p99 = zeros(T,1);
        for t = 1:T
            y = sort(deficits(t,:));
            k = max(1, ceil(0.99 * numel(y)));
            p99(t) = y(k);
        end
    end
    ypos = max(0, p99 - mu);  % upper-only whisker
    yneg = zeros(T,1);

    % threshold curve
    thresh = max(0, u_opt - u_min);

    % y-limit
    yMax = 1.15 * max([p99; thresh; eps]);

    figure('Color','w'); hold on; grid on; box on; ylim([0 yMax]);

    % violation zone
    xpoly = [tt; flipud(tt)];
    ypoly = [thresh; yMax*ones(T,1)];
    hshade = patch(xpoly, ypoly, [1 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.35);
    uistack(hshade,'bottom');

    % threshold line
    plot(tt, thresh, 'r--', 'LineWidth', 1.2);

    % bars (mean)
    bar(tt, mu, 'FaceColor',[0.3 0.6 0.9], 'EdgeColor','none');

    % upper-only 99th percentile whisker
    errorbar(tt, mu, yneg, ypos, 'k.', 'LineWidth', 1);

    xlabel('Time step');
    ylabel('Release Curtailment');
    title(titleStr);

    legend({'Violation Zone', ...
            'Violation Threshold', 'Mean Curtailment', '99th Percentile Spread'}, ...
            'Location','northwest');
end


function stats = analyzeFlowMinViolations(u_sim, u_min, unitLabel)
% analyzeFlowMinViolations
%
% Computes violation frequency, expected shortfall, conditional severity,
% 99th-percentile tail severity, etc., for simulated releases relative to u_min.
%
% Inputs:
%   u_sim     : T x nSim matrix (simulated curtailed turbine releases)
%   u_min     : scalar minimum release
%   unitLabel : optional string for printing titles (e.g. "Unit 1")
%
% Output (struct):
%   stats.p_violate_t   : T x 1 violation frequency per time step
%   stats.p_violate     : scalar overall violation frequency
%   stats.ES_t          : T x 1 expected shortfall per time step
%   stats.ES            : scalar overall expected shortfall
%   stats.ES_cond_t     : T x 1 conditional expected shortfall
%   stats.p99_short_t   : T x 1 99th-percentile shortfall
%
%   All metrics depend ONLY on u_sim and u_min — call AFTER simulations.

    if nargin < 3
        unitLabel = "Unit";
    end

    [T, nSim] = size(u_sim);

    % --- violation indicator + shortfall magnitude ---
    viol = (u_sim < u_min);                % T x nSim logical
    short = max(0, u_min - u_sim);         % T x nSim >= 0

    % ------------------------------------------------------------------
    % 1) Violation frequency
    % ------------------------------------------------------------------
    p_violate_t = mean(viol, 2);           % per time step
    p_violate   = mean(viol, 'all');       % overall

    % ------------------------------------------------------------------
    % 2) Expected shortfall (unconditional)
    % ------------------------------------------------------------------
    ES_t = mean(short, 2);                 % per time step
    ES   = mean(short, 'all');             % overall

    % ------------------------------------------------------------------
    % 3) Conditional severity (given a violation)
    % ------------------------------------------------------------------
    num_viol_t = sum(viol, 2);
    ES_cond_t  = zeros(T,1);
    nz = num_viol_t > 0;
    ES_cond_t(nz) = sum(short(nz,:), 2) ./ num_viol_t(nz);

    % ------------------------------------------------------------------
    % 4) Tail risk: 99th percentile shortfall
    % ------------------------------------------------------------------
    if exist('prctile', 'file')
        p99_short_t = prctile(short, 99, 2);
    else
        % fallback if prctile unavailable
        p99_short_t = zeros(T,1);
        for t = 1:T
            y = sort(short(t,:));
            k = max(1, ceil(0.99 * numel(y)));
            p99_short_t(t) = y(k);
        end
    end

    % ------------------------------------------------------------------
    % 5) Package into a nice struct
    % ------------------------------------------------------------------
    stats = struct();
    stats.p_violate_t = p_violate_t;
    stats.p_violate   = p_violate;
    stats.ES_t        = ES_t;
    stats.ES          = ES;
    stats.ES_cond_t   = ES_cond_t;
    stats.p99_short_t = p99_short_t;

    % ------------------------------------------------------------------
    % 6) Optional printed summary
    % ------------------------------------------------------------------
    fprintf('\nFlow-Min Violation Summary: %s\n', unitLabel);
    fprintf('Overall violation frequency: %.3f%%\n', 100*p_violate);
    fprintf('Overall expected shortfall (uncond.): %.4g\n', ES);
    fprintf('Mean conditional shortfall (given violate): %.4g\n', mean(ES_cond_t));
    fprintf('Max 99th-percentile shortfall: %.4g\n', max(p99_short_t));
end

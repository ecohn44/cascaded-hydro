function [V1_sim, V2_sim, u1_sim, u2_sim, p1_mean, p2_mean, deficits1, deficits2] = runBaseMonteCarloSims(sysparams, params, bounds, std_hat, X, savePath, printplot, curtail)
% runMonteCarloSims - Monte Carlo bound violation checker (FIXED POLICY, CLAMPED, WITH SIM HEAD)
%
%   Description:
%     (1) Added clamping to release/spill so volumes stay within bounds 
%     (2) Added generation recalculation from adjusted water release
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

    % Set random seed
    rng(0, 'twister');

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

        %  Extract inflow uncertainty (std) 
        std1 = std_hat(:,1);
        std2 = std_hat(:,2);

        for t = 1:T

            % Update standard deviation estimation for downstream units
            if t > 1 && curtail
                % Forecast conditional variance using GARCH-X
                var_hat = params.omega + 1e-6 + params.gamma*u1_sim(t-1, n);
                std2(t) = sqrt(var_hat);
            end 

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
            
            if curtail 
                % Reduce water release for volume boundaries
                [V1, u1, deficits1(t,n), violation_curtail_count1(t,n)] = curtail_deficit(V1_unc, u1_opt(t), Vmin1, Umin1);
                [V2, u2, deficits2(t,n), violation_curtail_count2(t,n)] = curtail_deficit(V2_unc, u2_opt(t), Vmin2, Umin2);
            else
                % Allow behavior outside volume constraints
                V1 = V1_unc;
                V2 = V2_unc;
            end


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

    % Function output
    p1_mean = mean(p1_sim, 2);
    p2_mean = mean(p2_sim, 2);

    % ===================================================================
    %  Report violation rates 
    % ===================================================================

    %  Violation audit (based on UNCLAMPED) 
    avg_per_timestep1 = mean(violation_count1, 2);
    avg_per_timestep2 = mean(violation_count2, 2);
    overall_avg1 = mean(avg_per_timestep1);
    overall_avg2 = mean(avg_per_timestep2);

    fprintf('\n  Monte Carlo Violation Report (UNCLAMPED audit) \n');
    fprintf('Reservoir 1: %.2f%%\n', 100*overall_avg1);
    fprintf('Reservoir 2: %.2f%%\n', 100*overall_avg2);

    % Calculate safety margin above Vmin
    volumeSafetyMargin(V1_sim, V2_sim, sysparams);

    %  Violation audit (based on CURTAILED)
    if curtail
        c_avg_per_timestep1 = mean(violation_curtail_count1, 2);
        c_avg_per_timestep2 = mean(violation_curtail_count2, 2);
        c_overall_avg1 = mean(c_avg_per_timestep1);
        c_overall_avg2 = mean(c_avg_per_timestep2);
    
        fprintf('\n  Monte Carlo Violation Report (CURTAILED audit) \n');
        fprintf('Reservoir 1: %.2f%%\n', 100*c_overall_avg1);
        fprintf('Reservoir 2: %.2f%%\n', 100*c_overall_avg2);
    
        reportFlowViolations(u1_sim, u2_sim, sysparams(1).min_ut, sysparams(2).min_ut);
        violationProbability(u1_sim, u2_sim, sysparams(1).min_ut, sysparams(2).min_ut);
        reportTotalRelease(u1_opt, u2_opt, u1_sim, u2_sim);
        reportTotalGeneration(p1_opt, p2_opt, p1_sim, p2_sim);
    end 

  
    % ===================================================================
    %  Visualizations
    % ===================================================================

    %{
    switch string(bounds)
        case "det",     bLabel = "Deterministic";
        case "icc",     bLabel = "Individual CC";
        case "jcc-bon", bLabel = "Bonferroni JCC";
        otherwise,      bLabel = char(bounds);
    end


    % Monte Carlo Volume Trajectories 

    plotReservoirMC(tt, V1_opt, Vmin1, Vmax1, ...
       sprintf('Reservoir 1 — %s', bLabel), 'Storage V_1 (p.u.)', 'reservoir1', savePath, bounds, printplot);
    plotReservoirMC(tt, V2_opt, Vmin2, Vmax2, ...
        sprintf('Reservoir 2 — %s', bLabel), 'Storage V_2 (p.u.)', 'reservoir2', savePath, bounds, printplot);

    %  Visualization of  Curtailment & Umin Violations

    plotCurtailment(deficits1, u1_opt, sysparams(1).min_ut, 'Unit 1 Curtailment');
    plotCurtailment(deficits2, u2_opt, sysparams(2).min_ut, 'Unit 2 Curtailment');

    %  Compare water release policies before and after curtailment 

    plotReleasePolicies(tt, u1_opt, u1_sim, sysparams(1).min_ut, sysparams(1).max_ut, deficits1, "Unit 1")
    plotReleasePolicies(tt, u2_opt, u2_sim, sysparams(2).min_ut, sysparams(2).max_ut, deficits2, "Unit 2")
    %}


end


function [V, u, take_turb, viol] = curtail_deficit(V_unc, u_opt, Vmin, u_min)

    % Record current volume and release
    V = V_unc; u = u_opt; take_turb = 0; viol = 0; 
    deficit_remain = 0;

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


function plotReservoirMC(tt, V_opt, Vmin, Vmax, titleStr, yLabelStr, ...
                         fileTag, savePath, bounds, printplot)

    % Mean trajectory 
    V_mean = mean(V_sim,2); 

    % IQR band
    pr1 = prctile(V_sim, [25 75], 2);   
    V_lower = pr1(:,1); 
    V_upper = pr1(:,2); 
    
    figure; hold on; grid on;

    % 1σ band
    fill([tt; flipud(tt)], [V_upper; flipud(V_lower)], ...
         [0.3 0.5 1.0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

    % Optimal and MC mean
    plot(tt, V_opt,  'r-',  'LineWidth', 2,   'DisplayName','Optimal');
    plot(tt, V_mean, 'b--', 'LineWidth', 1.5, 'DisplayName','MC mean');

    % Hard bounds
    yline(Vmin, '--r', 'Min bound');
    yline(Vmax, '--r', 'Max bound');

    xlabel('Hour');
    ylabel(yLabelStr);
    title(titleStr);

    legend('1\sigma band','Optimal','MC mean','Location','best');

    if printplot
        fname = sprintf('mc_%s_%s_fixed_clamped.png', fileTag, char(bounds));
        saveas(gcf, fullfile(savePath, fname));
    end
end


function plotReleasePolicies(tt, u_opt, u_mc, u_min, u_max, deficits, name)

    % Title
    titlestr = [name, ' Generation Release'];
    
    % Percentiles across simulations (row-wise)
    pr     = prctile(u_mc, [1 50 99], 2);
    u_p01  = pr(:,1);
    u_med  = pr(:,2);
    u_p99  = pr(:,3);
    
    % Mean curtailed profile (2nd bar)
    u_mean = mean(u_mc, 2);
    
    % --- Use deficits to detect if any curtailment happened at each t ---
    tol = 1e-9;                             % treat tiny values as zero
    curt_happened = any(deficits > tol, 2); % T x 1 logical
    
    % Zero whiskers where no curtailment occurred
    yneg = max(0, u_mean - u_p01);
    ypos = max(0, u_p99  - u_mean);
    yneg(~curt_happened) = 0;
    ypos(~curt_happened) = 0;
    
    % Grouped bars [Original, Curtailed]
    U_all = [u_opt, u_mean];
    
    figure('Color','w'); hold on; grid on; box on;
    b = bar(tt, U_all, 'grouped', 'BarWidth', 0.78);
    b(1).FaceColor = [0.2, 0.4, 0.8];   % soft dark blue
    b(2).FaceColor = [1.0, 1.0, 0.5];   % light yellow
    
    % Error bars (only where curtailment happened)
    x2 = b(2).XEndPoints;
    errorbar(x2, u_mean, yneg, ypos, 'k.', 'LineWidth', 1, 'CapSize', 3, ...
             'LineStyle','none');
    
    % Bounds (single legend entry)
    yline(u_min, '--r', 'LineWidth', 1.2, 'HandleVisibility','on');
    yline(u_max, '--r', 'LineWidth', 1.2, 'HandleVisibility','off');
    
    ylim([0, u_max*1.1]);
    xlabel('Time (h)'); ylabel('Generation Release');
    title(titlestr);
    set(gca, 'FontSize', 13);
    
    legend({'Original','Curtailed','Curtailed spread (p01–p99)','Bounds'}, ...
           'Location','southoutside','Orientation','horizontal','Box','off');
    
    % ===============================
    % Summary printout (using deficits mask)
    % For releases: worst case = LOWER tail → p01
    % ===============================
    total_opt  = sum(u_opt);
    total_mean = sum(u_mean);
    total_p01  = sum(u_p01);
    
    d_mean = total_mean - total_opt;
    d_p01  = total_p01  - total_opt;
    
    pct_mean = 100 * d_mean / max(total_opt, eps);
    pct_p01  = 100 * d_p01  / max(total_opt, eps);
    
    fprintf('\n Release Summary: %s\n', name);
    fprintf('Total (optimal policy):         %.4f\n', total_opt);
    fprintf('Total (curtailed, mean):        %.4f\n', total_mean);
    fprintf('   Δ (mean - opt):              %+8.4f  (%+6.2f%%)\n', d_mean, pct_mean);
    fprintf('Total (curtailed, 1st perc.):   %.4f\n', total_p01);
    fprintf('   Δ (p01 - opt):               %+8.4f  (%+6.2f%%)\n', d_p01,  pct_p01);

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


function reportFlowViolations(u1_sim, u2_sim, umin1, umin2)
% RELEASE_METRICS - Print Integrated Violation Index (IVI) summary
%                   for Unit 1, Unit 2, and the combined system.

    % Helper to compute IVI metrics for a single unit
    function [IVI_50, IVI_p01] = compute_IVI(u_sim, u_min)
        short_all = max(0, u_min - u_sim);  % violation magnitude per run
        A_all = sum(short_all, 1);          % total violation per run
        IVI_50 = median(A_all);             % median (typical run)
    
        u_p01 = prctile(u_sim, 1, 2);       % tail (1st percentile profile)
        short_p01 = max(0, u_min - u_p01);
        IVI_p01 = sum(short_p01);           % integrated tail violation
    end

    % Compute IVIs for Unit 1, Unit 2, and Combined
    [IVI50_1, IVIp01_1] = compute_IVI(u1_sim, umin1);
    [IVI50_2, IVIp01_2] = compute_IVI(u2_sim, umin2);

    % Total (weighted or direct sum across both units)
    IVI50_total  = IVI50_1 + IVI50_2;
    IVIp01_total = IVIp01_1 + IVIp01_2;

    fprintf('\n=================================================\n');
    fprintf(' Integrated Violation Index Summary (Flow-Min)\n');
    fprintf('=================================================\n');

    fprintf('\n  Unit 1:\n');
    fprintf('IVI (median run):        %.6f\n', IVI50_1);
    fprintf('IVI (1st percentile):    %.6f\n', IVIp01_1);

    fprintf('\n  Unit 2:\n');
    fprintf('IVI (median run):        %.6f\n', IVI50_2);
    fprintf('IVI (1st percentile):    %.6f\n', IVIp01_2);

    fprintf('\n  Total (Unit 1 + Unit 2):\n');
    fprintf('IVI (median run):        %.6f\n', IVI50_total);
    fprintf('IVI (1st percentile):    %.6f\n', IVIp01_total);

    fprintf('-------------------------------------------------\n');
end


function volumeSafetyMargin(V1_sim, V2_sim, sysparams)
% VOLUMESAFETYMARGIN  Computes and prints Volume Safety Margin Index (VSMI)

    % Extract parameters
    Vmin1 = sysparams(1).min_V;  Vmax1 = sysparams(1).max_V;
    Vmin2 = sysparams(2).min_V;  Vmax2 = sysparams(2).max_V;

    % Compute normalized distance above min
    safe1 = max(0, (V1_sim - Vmin1) ./ (Vmax1 - Vmin1));
    safe2 = max(0, (V2_sim - Vmin2) ./ (Vmax2 - Vmin2));

    % Mean safety margin across all time steps and simulations
    VSMI1 = mean(safe1(:));
    VSMI2 = mean(safe2(:));
    VSMI_total = (VSMI1 + VSMI2)/2;

    % Display results
    fprintf('\n=================================================\n');
    fprintf(' Volume Safety Margin Index (VSMI)\n');
    fprintf('=================================================\n');
    fprintf('  Unit 1: %.5f\n', VSMI1);
    fprintf('  Unit 2: %.5f\n', VSMI2);
    fprintf('  Total : %.5f (mean of units)\n', VSMI_total);
    fprintf('-------------------------------------------------\n');
end


function violationProbability(u1_sim, u2_sim, umin1, umin2)
% VIOLATIONPROBABILITY  Computes and prints Violation Probability Index (VPI)

    % Logical indicators of violation
    viol1 = (u1_sim < umin1);
    viol2 = (u2_sim < umin2);

    % Probability of violation across all runs/timesteps
    VPI1 = mean(viol1(:));
    VPI2 = mean(viol2(:));
    VPI_total = (VPI1 + VPI2)/2;

    % Display results
    fprintf('\n=================================================\n');
    fprintf(' Violation Probability Index (VPI)\n');
    fprintf('=================================================\n');
    fprintf('  Unit 1: %.4f%%\n', 100*VPI1);
    fprintf('  Unit 2: %.4f%%\n', 100*VPI2);
    fprintf('  Total : %.4f%% (mean of units)\n', 100*VPI_total);
    fprintf('-------------------------------------------------\n');
end


function reportTotalRelease(u1_opt, u2_opt, u1_sim, u2_sim)

    % Compute release profiles (CURTAILED) 
    u1_mean = mean(u1_sim, 2);
    u2_mean = mean(u2_sim, 2);

    % p01 (1st percentile) across simulations (pointwise)
    u1_p01 = prctile(u1_sim, 1, 2);
    u2_p01 = prctile(u2_sim, 1, 2);

    % Unit 1
    total_u1_opt  = sum(u1_opt);
    total_u1_mean = sum(u1_mean);
    total_u1_p01  = sum(u1_p01);

    denom1        = max(total_u1_opt, eps);
    change_u1_mean = 100 * (total_u1_mean - total_u1_opt) / denom1;
    change_u1_p01  = 100 * (total_u1_p01  - total_u1_opt) / denom1;

    % Unit 2
    total_u2_opt  = sum(u2_opt);
    total_u2_mean = sum(u2_mean);
    total_u2_p01  = sum(u2_p01);

    denom2        = max(total_u2_opt, eps);
    change_u2_mean = 100 * (total_u2_mean - total_u2_opt) / denom2;
    change_u2_p01  = 100 * (total_u2_p01  - total_u2_opt) / denom2;

    % System totals (Unit 1 + Unit 2)
    total_sys_opt  = total_u1_opt  + total_u2_opt;
    total_sys_mean = total_u1_mean + total_u2_mean;
    total_sys_p01  = total_u1_p01  + total_u2_p01;

    denom_sys       = max(total_sys_opt, eps);
    change_sys_mean = 100 * (total_sys_mean - total_sys_opt) / denom_sys;
    change_sys_p01  = 100 * (total_sys_p01  - total_sys_opt) / denom_sys;

    fprintf('\n=================================================\n');
    fprintf(' Release Comparison \n');
    fprintf('=================================================\n');
    fprintf(' Release Summary: System Total\n');
    fprintf('Total (optimal):         %.4f\n', total_sys_opt);
    fprintf('Total (mean curtailed):  %.4f  (Δ = %+8.4f, %+6.2f%%%%)\n', ...
            total_sys_mean, total_sys_mean - total_sys_opt, change_sys_mean);
    fprintf('Total (p01 curtailed):   %.4f  (Δ = %+8.4f, %+6.2f%%%%)\n', ...
            total_sys_p01,  total_sys_p01  - total_sys_opt, change_sys_p01);
    
    fprintf('\nRelease Summary: Unit 1\n');
    fprintf('Total (optimal):         %.4f\n', total_u1_opt);
    fprintf('Total (mean curtailed):  %.4f  (Δ = %+8.4f, %+6.2f%%%%)\n', ...
            total_u1_mean, total_u1_mean - total_u1_opt, change_u1_mean);
    fprintf('Total (p01 curtailed):   %.4f  (Δ = %+8.4f, %+6.2f%%%%)\n', ...
            total_u1_p01,  total_u1_p01  - total_u1_opt, change_u1_p01);
    
    fprintf('\nRelease Summary: Unit 2\n');
    fprintf('Total (optimal):         %.4f\n', total_u2_opt);
    fprintf('Total (mean curtailed):  %.4f  (Δ = %+8.4f, %+6.2f%%%%)\n', ...
            total_u2_mean, total_u2_mean - total_u2_opt, change_u2_mean);
    fprintf('Total (p01 curtailed):   %.4f  (Δ = %+8.4f, %+6.2f%%%%)\n', ...
            total_u2_p01,  total_u2_p01  - total_u2_opt, change_u2_p01);
    fprintf('-------------------------------------------------\n');
end


function reportTotalGeneration(p1_opt, p2_opt, p1_sim, p2_sim)

    % Curtailed means (per time step)
    p1_mean = mean(p1_sim, 2);
    p2_mean = mean(p2_sim, 2);

    % Lower-tail (worst) generation per time step: p01
    p1_p01  = prctile(p1_sim, 1, 2);
    p2_p01  = prctile(p2_sim, 1, 2);

    % Totals (sum over time) 
    tot_p1_opt  = sum(p1_opt);
    tot_p2_opt  = sum(p2_opt);
    tot_p1_mean = sum(p1_mean);
    tot_p2_mean = sum(p2_mean);
    tot_p1_p01  = sum(p1_p01);
    tot_p2_p01  = sum(p2_p01);

    % Per-unit deltas and percents 
    d1_mean   = tot_p1_mean - tot_p1_opt;
    d1_p01    = tot_p1_p01  - tot_p1_opt;
    pct1_mean = 100 * d1_mean / max(tot_p1_opt, eps);
    pct1_p01  = 100 * d1_p01 / max(tot_p1_opt, eps);

    d2_mean   = tot_p2_mean - tot_p2_opt;
    d2_p01    = tot_p2_p01  - tot_p2_opt;
    pct2_mean = 100 * d2_mean / max(tot_p2_opt, eps);
    pct2_p01  = 100 * d2_p01 / max(tot_p2_opt, eps);

    % Combined totals (system) 
    total_p_opt  = tot_p1_opt  + tot_p2_opt;
    total_p_mean = tot_p1_mean + tot_p2_mean;
    total_p_p01  = tot_p1_p01  + tot_p2_p01;

    d_mean_tot   = total_p_mean - total_p_opt;
    d_p01_tot    = total_p_p01  - total_p_opt;

    pct_mean_tot = 100 * d_mean_tot / max(total_p_opt, eps);
    pct_p01_tot  = 100 * d_p01_tot  / max(total_p_opt, eps);

    % Total system summary
    fprintf('\n=================================================\n');
    fprintf(' Total Generation Comparison\n');
    fprintf('=================================================\n');
    fprintf('Total (optimal):                %.4f\n', total_p_opt);
    fprintf('Total (mean curtailed):         %.4f  (Δ = %+8.4f, %+6.2f%%%%)\n', ...
            total_p_mean, d_mean_tot, pct_mean_tot);
    fprintf('Total (p01 curtailed):          %.4f  (Δ = %+8.4f, %+6.2f%%%%)\n', ...
            total_p_p01,  d_p01_tot,  pct_p01_tot);

    % Per-unit summaries
    fprintf('\n Generation Summary: Unit 1\n');
    fprintf('Total (optimal policy):         %.4f\n', tot_p1_opt);
    fprintf('Total (curtailed, mean):        %.4f\n', tot_p1_mean);
    fprintf('   Δ (mean - opt):              %+8.4f  (%+6.2f%%%%)\n', d1_mean,  pct1_mean);
    fprintf('Total (curtailed, p01):         %.4f\n', tot_p1_p01);
    fprintf('   Δ (p01 - opt):               %+8.4f  (%+6.2f%%%%)\n', d1_p01,   pct1_p01);

    fprintf('\n Generation Summary: Unit 2\n');
    fprintf('Total (optimal policy):         %.4f\n', tot_p2_opt);
    fprintf('Total (curtailed, mean):        %.4f\n', tot_p2_mean);
    fprintf('   Δ (mean - opt):              %+8.4f  (%+6.2f%%%%)\n', d2_mean,  pct2_mean);
    fprintf('Total (curtailed, p01):         %.4f\n', tot_p2_p01);
    fprintf('   Δ (p01 - opt):               %+8.4f  (%+6.2f%%%%)\n', d2_p01,   pct2_p01);
    fprintf('-------------------------------------------------\n');
end
function [V1, V2] = runMonteCarloSims(sysparams, bounds, std_hat, X, savePath, printplot)
% runMonteCarloSims - Monte Carlo bound violation checker (FIXED POLICY)
%
%   Always uses the policy exactly as solved: u_t = u_opt(t), s_t = s_opt(t).
%   Proper mass balance: V_t = V_{t-1} + q_t - u_t - s_t
%
%   Inputs:
%       sysparams : struct array (size 2) with fields .min_V, .max_V, .V0
%       bounds    : type of volume bounds used (for labeling only)
%       std_hat   : Tx2 matrix of forecast inflow std deviations (same UNITS as q)
%       X         : Tx10 matrix of deterministic simulation results
%                   Columns:
%                     1=V1, 2=p1, 3=u1, 4=s1, 5=q1,
%                     6=V2, 7=p2, 8=u2, 9=s2, 10=q2
%       savePath  : folder for plots
%       printplot : true/false to save figures
%
%   Outputs:
%       V1, V2    : TxN simulated volume trajectories (for reference)

    % Monte Carlo samples
    nSim = 1000;

    % Extract T
    T = size(X,1);

    % Volume bounds and initials
    Vmin1 = sysparams(1).min_V; Vmax1 = sysparams(1).max_V; V0_1 = sysparams(1).V0;
    Vmin2 = sysparams(2).min_V; Vmax2 = sysparams(2).max_V; V0_2 = sysparams(2).V0;

    % Optimal (solved) trajectories / policy
    V1_opt = X(:,1);  u1_opt = X(:,3);  s1_opt = X(:,4);  q1_mean = X(:,5);
    V2_opt = X(:,6);  u2_opt = X(:,8);  s2_opt = X(:,9);  q2_mean = X(:,10);

    % Std devs (MUST be in same units as q_mean)
    std1 = std_hat(:,1);
    std2 = std_hat(:,2);

    % Storage for MC runs
    V1 = zeros(T, nSim);
    V2 = zeros(T, nSim);
    violation_count1 = false(T, nSim);
    violation_count2 = false(T, nSim);

    % --- Monte Carlo propagation (FIXED POLICY) ---
    for n = 1:nSim
        V1_prev = V0_1; V2_prev = V0_2;

        for t = 1:T
            % Sample inflows
            inflow1 = q1_mean(t) + std1(t)*randn;
            inflow2 = q2_mean(t) + std2(t)*randn;

            % Use the solved policy exactly
            u1 = u1_opt(t); s1 = s1_opt(t);
            u2 = u2_opt(t); s2 = s2_opt(t);

            % Proper mass balance **includes spill**
            V1_t = V1_prev + inflow1 - u1 - s1;
            V2_t = V2_prev + inflow2 - u2 - s2;

            % Record
            V1(t,n) = V1_t; V2(t,n) = V2_t;

            % Violations
            violation_count1(t,n) = (V1_t < Vmin1) || (V1_t > Vmax1);
            violation_count2(t,n) = (V2_t < Vmin2) || (V2_t > Vmax2);

            % Advance
            V1_prev = V1_t; V2_prev = V2_t;
        end
    end

    % Violation rates
    avg_per_timestep1 = mean(violation_count1, 2);
    overall_avg1 = mean(avg_per_timestep1);
    avg_per_timestep2 = mean(violation_count2, 2);
    overall_avg2 = mean(avg_per_timestep2);

    fprintf('\n ---Monte Carlo Violation Report---\n');
    fprintf('Reservoir 1: %.2f%%\n', 100*overall_avg1);
    fprintf('Reservoir 2: %.2f%%\n', 100*overall_avg2);

    % MC envelopes for visualization (1σ band)
    V1_mean = mean(V1, 2);  V1_std = std(V1, 0, 2);
    V2_mean = mean(V2, 2);  V2_std = std(V2, 0, 2);

    V1_upper = V1_mean + V1_std; V1_lower = V1_mean - V1_std;
    V2_upper = V2_mean + V2_std; V2_lower = V2_mean - V2_std;

    % Pretty label for bounds
    switch bounds
        case "det",     bLabel = "Deterministic";
        case "icc",     bLabel = "Individual CC";
        case "jcc-bon", bLabel = "Bonferroni JCC";
        otherwise,      bLabel = char(bounds);
    end

    % --- Plot Reservoir 1 ---
    figure; hold on; grid on;
    fill([1:T, fliplr(1:T)], [V1_upper', fliplr(V1_lower')], ...
         [0.3 0.5 1.0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(V1_opt, 'r-', 'LineWidth', 2, 'DisplayName','Optimal V1');
    plot(V1_mean, 'b--', 'LineWidth', 1.5, 'DisplayName','MC mean');
    yline(Vmin1, '--r', 'Min bound'); yline(Vmax1, '--r', 'Max bound');
    xlabel('Hour'); ylabel('Storage V_1 (m^3)');
    title(sprintf('Reservoir 1 — %s', bLabel));
    legend('1\sigma band','Optimal','MC mean','Location','best');
    if printplot, saveas(gcf, fullfile(savePath, ['mc_reservoir1_' char(bounds) '_fixed.png'])); end

    % --- Plot Reservoir 2 ---
    figure; hold on; grid on;
    fill([1:T, fliplr(1:T)], [V2_upper', fliplr(V2_lower')], ...
         [0.3 0.5 1.0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(V2_opt, 'r-', 'LineWidth', 2, 'DisplayName','Optimal V2');
    plot(V2_mean, 'b--', 'LineWidth', 1.5, 'DisplayName','MC mean');
    yline(Vmin2, '--r', 'Min bound'); yline(Vmax2, '--r', 'Max bound');
    xlabel('Hour'); ylabel('Storage V_2 (m^3)');
    title(sprintf('Reservoir 2 — %s', bLabel));
    legend('1\sigma band','Optimal','MC mean','Location','best');
    if printplot, saveas(gcf, fullfile(savePath, ['mc_reservoir2_' char(bounds) '_fixed.png'])); end
end

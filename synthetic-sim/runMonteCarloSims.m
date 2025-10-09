function [V1, V2] = runMonteCarloSims(sysparams, bounds, std_hat, X)
% runMonteCarloCheck - Monte Carlo bound violation checker for reservoir simulation
%
%   Inputs:
%       sysparams : struct array (size 2) with fields .min_V and .max_V
%       std_hat   : Tx2 matrix of forecast inflow standard deviations
%       X         : Tx10 matrix of deterministic simulation results
%                   Columns: 
%                     1=V1, 2=p1, 3=u1, 4=s1, 5=q1, 
%                     6=V2, 7=p2, 8=u2, 9=s2, 10=q2

    nSim = 5000;                % Monte Carlo samples
    T = size(X,1);
    tvec = 1:T;

    % Extract bounds
    Vmin1 = sysparams(1).min_V;
    Vmax1 = sysparams(1).max_V;
    Vmin2 = sysparams(2).min_V;
    Vmax2 = sysparams(2).max_V;

    % Optimial Trajectories
    V1_opt = X(:,1);  u1_opt = X(:,3);  q1_mean = X(:,5);
    V2_opt = X(:,6);  u2_opt = X(:,8);  q2_mean = X(:,10);

    std1 = std_hat(:,1);
    std2 = std_hat(:,2);

    % Initialize
    V1 = zeros(T, nSim);
    V2 = zeros(T, nSim);
    V1(1,:) = V1_opt(1);
    V2(1,:) = V2_opt(1);

    % Monte Carlo propagation
    for t = 2:T
        inflow1 = q1_mean(t) + randn(1, nSim) * std1(t);
        inflow2 = q2_mean(t) + randn(1, nSim) * std2(t);

        V1(t,:) = V1(t-1,:) + inflow1 - u1_opt(t);
        V2(t,:) = V2(t-1,:) + inflow2 - u2_opt(t);
    end

    % Violation rates
    viol_low1  = mean(V1(:) < Vmin1);
    viol_high1 = mean(V1(:) > Vmax1);
    viol_low2  = mean(V2(:) < Vmin2);
    viol_high2 = mean(V2(:) > Vmax2);

    % Display
    fprintf('\n=== Monte Carlo Violation Check ===\n');
    fprintf('Reservoir 1: below min = %.2f%%, above max = %.2f%%\n', 100*viol_low1, 100*viol_high1);
    fprintf('Reservoir 2: below min = %.2f%%, above max = %.2f%%\n', 100*viol_low2, 100*viol_high2);

    % Compute ensemble stats
    meanV1 = mean(V1, 2);
    stdV1  = std(V1, 0, 2);
    meanV2 = mean(V2, 2);
    stdV2  = std(V2, 0, 2);

    switch bounds
        case "det"
            bLabel = "Deterministic";
        case "icc"
            bLabel = "Individual CC";
        case "jcc-bon"
            bLabel = "Bonferroni JCC";
    end

    % === Plot Reservoir 1 ===
    figure;
    hold on; grid on;
    fill([tvec fliplr(tvec)], [meanV1-stdV1; flipud(meanV1+stdV1)]', [0.7 0.8 1], ...
         'FaceAlpha', 0.4, 'EdgeColor', 'none');
    plot(tvec, V1_opt, 'b-', 'LineWidth', 2);
    yline(Vmin1, '--r', 'Min bound');
    yline(Vmax1, '--r', 'Max bound');
    xlabel('Time step');
    ylabel('Reservoir 1 Volume');
    title('Reservoir 1: Storage trajectories (mean ±1σ)');
    legend('MC range (±1σ)', bLabel, 'Bounds');

    % === Plot Reservoir 2 ===
    figure;
    hold on; grid on;
    fill([tvec fliplr(tvec)], [meanV2-stdV2; flipud(meanV2+stdV2)]', [0.7 0.8 1], ...
         'FaceAlpha', 0.4, 'EdgeColor', 'none');
    plot(tvec, V2_opt, 'b-', 'LineWidth', 2);
    yline(Vmin2, '--r', 'Min bound');
    yline(Vmax2, '--r', 'Max bound');
    xlabel('Time step');
    ylabel('Reservoir 2 Volume');
    title('Reservoir 2: Storage trajectories (mean ±1σ)');
    legend('MC range (±1σ)', bLabel, 'Bounds');
end

function [V1, V2] = runMonteCarloSims(sysparams, bounds, std_hat, X, savePath, printplot)
% runMonteCarloCheck - Monte Carlo bound violation checker for reservoir simulation
%
%   Inputs:
%       sysparams : struct array (size 2) with fields .min_V and .max_V
%       bounds    : type of volume bounds used
%       std_hat   : Tx2 matrix of forecast inflow standard deviations
%       X         : Tx10 matrix of deterministic simulation results
%                   Columns: 
%                     1=V1, 2=p1, 3=u1, 4=s1, 5=q1, 
%                     6=V2, 7=p2, 8=u2, 9=s2, 10=q2


    % Monte Carlo samples
    nSim = 1000;  
    % Extract T from simulations 
    T = size(X,1);

    % Extract Volume bounds
    Vmin1 = sysparams(1).min_V;
    Vmax1 = sysparams(1).max_V;
    Vmin2 = sysparams(2).min_V;
    Vmax2 = sysparams(2).max_V;

    % Optimial Trajectories
    V1_opt = X(:,1);  u1_opt = X(:,3);  q1_mean = X(:,5);
    V2_opt = X(:,6);  u2_opt = X(:,8);  q2_mean = X(:,10);

    % Extract estimated std deviations 
    std1 = std_hat(:,1);
    std2 = std_hat(:,2);

    % Initialize storage vectors
    V1 = zeros(T, nSim);
    V2 = zeros(T, nSim);
    violation_count1 = zeros(T, nSim);
    violation_count2 = zeros(T, nSim);

    % Monte Carlo propagation
    for n = 1:nSim
        for t = 1:T
            % Sample inflow
            inflow1 = normrnd(q1_mean(t), std1(t));
            inflow2 = normrnd(q2_mean(t), std2(t));
            
            % Calculuate volume trajectory
            if t == 1
                V1(t,n) = sysparams(1).V0 + inflow1 - u1_opt(t);
                V2(t,n) = sysparams(2).V0 + inflow2 - u2_opt(t);
            else
                V1(t,n) = V1(t-1,n) + inflow1 - u1_opt(t);
                V2(t,n) = V2(t-1,n) + inflow2 - u2_opt(t);
            end 

            % Violation count per time step (both bounds at once)
            if V1(t,n) < Vmin1 || V1(t,n) > Vmax1
                % Raise violation count 
                violation_count1(t,n) = 1;
            end

             % Violation count per time step (both bounds at once)
            if V2(t,n) < Vmin2 || V2(t,n) > Vmax2
                % Raise violation count 
                violation_count2(t,n) = 1;
            end
        end
    end

    % Calculate critical area (where v(t) ~ Vmin)
    % critical_timesteps1 = V1_opt <= min(V1_opt);
    % critical_timesteps2 = V2_opt <= min(V2_opt);

    % Average violations in the critical area
    avg_per_timestep1 = mean(violation_count1, 2);
    overall_avg1 = mean(avg_per_timestep1);
    avg_per_timestep2 = mean(violation_count2, 2);
    overall_avg2 = mean(avg_per_timestep2);

    % Display percent violations
    fprintf('\n ---Monte Carlo Violation Report---\n');
    fprintf('Reservoir 1: %.2f%%\n', 100*overall_avg1);
    fprintf('Reservoir 2: %.2f%%\n', 100*overall_avg2);

    % Calculate mean and std across simulations
    V1_mean = mean(V1, 2);  % Mean across simulations for each timestep
    V1_std = std(V1, 0, 2);  % Std deviation across simulations
    V2_mean = mean(V2, 2);  
    V2_std = std(V2, 0, 2);  
    
    % Upper and lower bounds (1 standard deviation)
    V1_upper = V1_mean + V1_std;
    V1_lower = V1_mean - V1_std;
    V2_upper = V2_mean + V2_std;
    V2_lower = V2_mean - V2_std;

    switch bounds
        case "det"
            bLabel = "Deterministic";
        case "icc"
            bLabel = "Individual CC";
        case "jcc-bon"
            bLabel = "Bonferroni JCC";
    end

    % Plot Reservoir 1 Monte Carlo Sims

    figure;
    hold on; grid on;
    fill([1:T, fliplr(1:T)], [V1_upper', fliplr(V1_lower')], ...
         'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on;
    plot(V1_opt, 'r-', 'LineWidth', 2);
    plot(V1_mean, 'b--', 'LineWidth', 1.5);
    yline(Vmin1, '--r', 'Min bound');
    yline(Vmax1, '--r', 'Max bound');
    xlabel('Time Step');
    ylabel('Storage Volume');
    legend('1 Std Dev', 'V1_{opt}', 'MC Mean', 'Location', 'best');
    title('Optimal Trajectory vs Monte Carlo Simulations');
    grid on;
    if printplot
        saveas(gcf, fullfile(savePath, ['mc_reservoir1_' char(bounds) '.png']));
    end

    % Plot Reservoir 2 Monte Carlo Sims
    figure;
    hold on; grid on;
    fill([1:T, fliplr(1:T)], [V2_upper', fliplr(V2_lower')], ...
         'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on;
    plot(V2_opt, 'r-', 'LineWidth', 2);
    plot(V2_mean, 'b--', 'LineWidth', 1.5);
    yline(Vmin2, '--r', 'Min bound');
    yline(Vmax2, '--r', 'Max bound');
    xlabel('Time Step');
    ylabel('Storage Volume');
    legend('1 Std Dev', 'V2_{opt}', 'MC Mean', 'Location', 'best');
    title('Optimal Trajectory vs Monte Carlo Simulations');
    grid on;
    if printplot
        saveas(gcf, fullfile(savePath, ['mc_reservoir2_' char(bounds) '.png']));
    end

end

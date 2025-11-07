%% Plot Driver for DIU vs DIU Comparison
% ========================================================================
% Author: Eliza Cohn
% Description: Overlays policy behavior derived during optimization sims
% Inputs: For uncertainty framework
%       X: Trajectory of optimal deicision variables 
%       std_hat: Standard deviation of inflow 
%       V_eff: Effective volume boundaries 
%       sysparams: System parameters
% ========================================================================

clc; clear; close all;

% Load optimal trajectories for each policy (for units n > 1)
M1 = load('./results/results_unit2_det.mat');
M2 = load('./results/results_unit2_diu.mat');
M3 = load('./results/results_unit2_ddu.mat');

% Plot settings
season = M1.season;
s = M1.sysparams;
opts = struct('titleTag',season);
T  = M1.T;
tt = (1:T)';
Vmin  = s(2).min_V;
Vmax  = s(2).max_V;

savePath = './results/';
if ~exist(savePath,'dir')
    mkdir(savePath);
end


% =====================================================
% FIGURE 1: Standard Deviation Comparison for M2 and M3
% =====================================================

f1 = figure('Position',[100 100 850 350]); hold on; grid on;
plot(tt, M2.std_hat(:,2), 'LineWidth',3, 'DisplayName','M2 \sigma_t');
plot(tt, M3.std_hat(:,2), 'LineWidth',3, 'DisplayName','M3 \sigma_t');
ylabel('Standard Deviation');
xlabel('Time (h)')
title('Unit 02 Inflow Uncertainty');
legend('Location','northeast');

saveas(f1, fullfile(savePath, 'u2_sigma_diu_vs_ddu.png'));

% =====================================================
% FIGURE 2: Effective Volume Boundaries
% =====================================================
% V_eff format is: [V1_max V1_min V2_max V2_min]

% Extract Volume Bounds/Trajectories for Unit 2
V2max_M1 = M1.V_eff(:,3);
V2min_M1 = M1.V_eff(:,4);
V2_M1 = M1.X(:,6);

V2max_M2 = M2.V_eff(:,3);
V2min_M2 = M2.V_eff(:,4);
V2_M2 = M2.X(:,6);

V2max_M3 = M3.V_eff(:,3);
V2min_M3 = M3.V_eff(:,4);
V2_M3 = M3.X(:,6);

% Base MATLAB default colors
blueBound  = [0 0 1];         % blue
greenBound = [0 1 0];         % green
redBound   = [1 0 0];         % red

% Trajectory Colors
blueTraj  = 0.6 * blueBound;
greenTraj = 0.6 * greenBound;
redTraj   = 0.6 * redBound;

f2 = figure('Position',[100 480 1100 480]); hold on; grid on;

% M1: Deterministic bounds and policy
plot(tt, V2min_M1, 'LineWidth',3, 'LineStyle','-.', 'Color',blueBound, 'DisplayName','M1^{Bounds}');
plot(tt, V2max_M1, 'LineWidth',3, 'LineStyle','-.', 'Color',blueBound, 'HandleVisibility','off');
plot(tt, V2_M1,    'LineWidth',2, 'Color',blueTraj, 'DisplayName','M1^{Trajectory}');

% M2: DIU bounds and policy
plot(tt, V2min_M2, 'LineWidth',3, 'LineStyle','-.', 'Color',greenBound, 'DisplayName','M2^{Bounds}');
plot(tt, V2max_M2, 'LineWidth',3, 'LineStyle','-.', 'Color',greenBound, 'HandleVisibility','off');
plot(tt, V2_M2,    'LineWidth',2, 'Color',greenTraj, 'DisplayName','M2^{Trajectory}');

% M3: DDU bounds and policy
plot(tt, V2min_M3, 'LineWidth',3, 'LineStyle','-.', 'Color',redBound, 'DisplayName','M3^{Bounds}');
plot(tt, V2max_M3, 'LineWidth',3, 'LineStyle','-.', 'Color',redBound, 'HandleVisibility','off');
plot(tt, V2_M3,    'LineWidth',2, 'Color',redTraj, 'DisplayName','M3^{Trajectory}');

% ylim([0,0.06]);
ylabel('SoC Bounds');
xlabel('Time (h)')
set(gca, 'FontSize', 14);        % axes labels, ticks
set(findall(gcf,'Type','text'), 'FontSize', 14);  % all text objects
legend('Location','best', 'FontSize', 13);
title('SoC Trajectory and Effective Bounds', 'FontSize', 16);

saveas(f2, fullfile(savePath, 'u2_soc_diu_vs_ddu.png'));

% =====================================================
% FIGURE 3: Power Production (Normalized)
% =====================================================

% Extract power dispatch for each policy
p2_M1 = M1.X(:,7);
p2_M2 = M2.X(:,7);
p2_M3 = M3.X(:,7);

P_all = [p2_M1, p2_M2, p2_M3];

f3 = figure('Position',[100 480 1100 480]); hold on; grid on;

% Plot power dispatch 
b = bar(tt, P_all, 'grouped');

% Apply MATLAB default colors
b(1).FaceColor = [0 0 1];    % blue
b(2).FaceColor = [0 1 0];    % green
b(3).FaceColor = [1 0 0];    % red

ylabel('Normalized Power (p.u.)');
xlabel('Time (h)')
set(gca, 'FontSize', 14);        % axes labels, ticks
set(findall(gcf,'Type','text'), 'FontSize', 14);  % all text objects
legend({'M1', 'M2', 'M3'}, 'Location', 'best', 'FontSize', 13);
title('Unit 02 Normalized Power Dispatch', 'FontSize', 16);

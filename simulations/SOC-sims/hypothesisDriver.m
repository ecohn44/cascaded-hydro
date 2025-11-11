%% Plot Driver for DIU vs DIU Comparison
% ========================================================================
% Author: Eliza Cohn
% Description: Test perforamnce of various policies under uncertainty
% frameworks in a Monte Carlo simulation setting 
% 
% Inputs: For uncertainty framework
%       X: Trajectory of optimal deicision variables 
%       std_hat: Standard deviation of inflow 
%       V_eff: Effective volume boundaries 
%       sysparams: System parameters
% ========================================================================

clc; clear; close all;

% Add shared functions to file path 
thisFilePath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisFilePath, '..', 'functions')));

% Load optimal trajectories for each policy (for units n > 1)
M1 = load('./results/results_unit2_det.mat');
M2 = load('./results/results_unit2_diu.mat');
M3 = load('./results/results_unit2_ddu.mat');
T  = M1.T;
tt = (1:T)';

% For Hypothesis 02, load in std_hat from DDU run 
std_hat = M3.std_hat; 

% Load simulation settings 
simSettings.bounds = "jcc-bon";
s = M1.sysparams;
printplot = false; 
path = ""; % Don't save pngs for now

% Run MC Sims for M1 under DDU assumption 
%[M1_V1, M1_V2, M1_u1, M2_u2, M1_p1_mean, M1_p2_mean, d1, d2] = runMonteCarloSims(s, simSettings.bounds, std_hat, M1.X, path, printplot);

% Run MC Sims for M2 under DDU assumption 
% [M2_V1, M2_V2, M2_p1_mean, M2_p2_mean] = runMonteCarloSims(s, simSettings.bounds, std_hat, M2.X, path, printplot);

% Run MC Sims for M3 under DDU assumption 
[M3_V1, M3_V2, M3_p1_mean, M3_p2_mean] = runMonteCarloSims(s, simSettings.bounds, std_hat, M3.X, path, printplot);


%{
% =====================================================
% FIGURE 1: Clamped Power Production (Normalized)
% =====================================================

% Extract power dispatch for each policy
P_all = [M1_p2_mean, M2_p2_mean, M3_p2_mean];

f1 = figure('Position',[100 480 1100 480]); hold on; grid on;

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

% =====================================================
% FIGURE 2: Curtailed Power Production vs Optimal Power Production 
% =====================================================

P_M1_all = [M1.X(:,7), M1_p2_mean];
f2 = figure('Position',[100 480 1100 480]); hold on; grid on;

% Plot power dispatch 
b = bar(tt, P_M1_all, 'grouped');

ylabel('Normalized Power (p.u.)');
xlabel('Time (h)')
set(gca, 'FontSize', 14);        % axes labels, ticks
set(findall(gcf,'Type','text'), 'FontSize', 14);  % all text objects
legend({'M1 Original', 'M1 Clamped'}, 'Location', 'best', 'FontSize', 13);
title('Unit 02 Normalized Power Dispatch', 'FontSize', 16);

% =====================================================
% FIGURE 3: Curtailed Power Production vs Optimal Power Production 
% =====================================================

%Unit 1
P_M1_all = [M1.X(:,2), M1_p1_mean];
f3 = figure('Position',[100 480 1100 480]); hold on; grid on;

% Plot power dispatch 
b = bar(tt, P_M1_all, 'grouped');

ylabel('Normalized Power (p.u.)');
xlabel('Time (h)')
set(gca, 'FontSize', 14);        % axes labels, ticks
set(findall(gcf,'Type','text'), 'FontSize', 14);  % all text objects
legend({'M1 Original', 'M1 Clamped'}, 'Location', 'best', 'FontSize', 13);
title('Unit 01 Normalized Power Dispatch', 'FontSize', 16);
%}
% clear; clc; close all;

% --- Load results from the universal results folder ---
S_diu = load('./results/results_unit2_diu.mat');
S_ddu = load('./results/results_unit2_ddu.mat');

% --- Plot options ---
opts = struct('showBands', true, 'showDelta', false, 'titleTag', 'dry');

% --- Call the overlay plotting function ---
plotDIUvDDU( ...
    S_diu.X, S_diu.std_hat, ...
    S_ddu.X, S_ddu.std_hat, ...
    S_ddu.sysparams, './plots', opts, ...
    S_diu.U_eff, S_ddu.U_eff);

%% Author: Eliza Cohn
% Date: November 2025
% Description: Main driver for cascaded hydropower simulations 
% Paper: Optimization of Cascaded Hydroelectric Systems under DDU

clear; clc; close all;
addpath('/Library/gurobi1202/macos_universal2/matlab');

% Add shared functions to file path 
thisFilePath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisFilePath, '..', 'functions')));

%% ========================================================================
% SECTION 1: DATA LOADING AND PARAMETER DEFINITION
% ========================================================================

% Static Parameters (PowerProd)
eta = .9;          % efficiency of release-energy conversion
rho_w = 1000;      % density of water [kg/m^3]
g = 9.8;           % acceleration due to gravity [m/s^2]
% c = eta*rho_w*g/3.6e9; % power prod coefficient
c = 1;
N = 20;             % number of sub-intervals for piecewise linear approx
n = 4;              % number of units in cascaded network 

% Load inflow data 
[modelparams, sysparams, droughtparams] = dataload(n, N);

%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings (season, drought type, linear approximation, uncertainty, bounds)
simSettings = initSimSettings("dry", "extended", "pwl", "diu", "jcc-bon");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Extract drought simulation mode
droughtparams = droughtparams(strcmp({droughtparams.mode}, simSettings.drought));

% Date range settings 
D = 7;                        % Simulation duration in days
T = 24*D;                     % Number of simulation hours
lag = 1;                      % Number of lag terms in OLS model

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: STREAMFLOW BEHAVIOR
% ========================================================================

% Matrix for local streamflow 
q = zeros(T+lag, n);

% Base inflow
q0 = 0.075;   

% Simulate drought event 
q(:,1) = droughtSimulator(q0, T, lag, simSettings.season, droughtparams.mode, droughtparams);

%% ========================================================================
% SECTION 4: OPTIMIZATION FRAMEWORK
% ========================================================================

scale = 1; % Scale safety bounds 
% scale_ddu = 25; % Scale gamma 

[model, obj, X, std_hat, V_eff] = genOptimization(T, N, c, q, lag, scale, ...
    simSettings.framework, simSettings.bounds, modelparams, sysparams);


%% ========================================================================
% SECTION 5: PLOTTING
% ========================================================================

% Toggle for creating folder and plotting
make_dir = false;
printplot = false; 

% Make plot directory for current simulation run 
if make_dir
    dir_path = "./plots/";
    stamp = datestr(now,'mm-dd-yyyy HH.MM.SS');
    path = fullfile(dir_path, stamp + " " + simSettings.season + " " ...
        + simSettings.framework + ...
        " T=" + string(T));
    mkdir(path)
end

% Plot simulation behavior for all units
simPlots(path, X, sysparams, T, c, printplot);

% Store simulation results 
results_dir = "./results/";
for i = 1:numel(sysparams)
    sp = sysparams(i);
    fname = sprintf('results_unit%d_%s.mat', ...
        sp.unit, lower(simSettings.framework));

    season = simSettings.season;  

    save(fullfile(results_dir, fname), ...
        'X', 'V_eff', 'std_hat', 'q', 'sysparams', ...
        'T', 'c', 'lag', 'season', '-v7');
end


fprintf('Simulation complete.\n');


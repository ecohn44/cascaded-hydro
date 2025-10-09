%% Author: Eliza Cohn
% Date: October 2025
% Description: Main driver for cascaded hydropower simulations 
% Paper: Optimization of Cascaded Hydroelectric Systems under DDU

clear; clc; close all;
addpath('/Library/gurobi1202/macos_universal2/matlab');

%% ========================================================================
% SECTION 1: DATA LOADING AND PARAMETER DEFINITION
% ========================================================================

% Static Parameters (PowerProd)
eta = .9;          % efficiency of release-energy conversion
rho_w = 1000;      % density of water [kg/m^3]
g = 9.8;           % acceleration due to gravity [m/s^2]
c = eta*rho_w*g/3.6e9; % power prod coefficient
N = 20;             % number of sub-intervals for piecewise linear approx

% Load inflow data
[inflow, modelparams, sysparams] = dataload(N);

%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings (season, linear approximation, uncertainty, bounds)
simSettings = initSimSettings("dry", "pwl", "diu", "jcc-bon");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));
modelparams.rho = 0.1135; % Calculated offline between (q1_hist, s1 + u1) 

% Date range settings 
D = 5;                        % Simulation duration in days
% T = 12 + 24*D;                % Number of simulation hours
T = 60;
lag = 1;                      % Number of lag terms in OLS model
year = 2022;                  % Simulation year

% Compute simulation daterange and inflow series
sim_center_date = datetime(year, 1, 1) + days(modelparams.center_day - 1);
start_date = sim_center_date - hours(T/2) - hours(lag);
end_date   = sim_center_date + hours(T/2 - 1);
inflow_s = inflow(inflow.datetime >= start_date & inflow.datetime <= end_date, :);

% Extract historic inflow timeseries [m3/hr]
q = [inflow_s.bon_inflow_m3hr];

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% ========================================================================

[model, obj, X, std_hat] = optimization(T, N, c, q, lag, ...
    simSettings.framework, simSettings.bounds, modelparams, sysparams);

% Extract q2 reference inflow
q(:,2) = [0; X(:,3) + X(:,4)];

%% ========================================================================
% SECTION 4: PLOTTING
% ========================================================================

% Toggle for creating folder and plotting
make_dir = true; % Set to true to enable directory creation and plotting
printplot = false; 

if make_dir
    dir_path = "./plots/";
    stamp = datestr(now,'mm-dd-yyyy HH.MM.SS');
    path = fullfile(dir_path, stamp + " " + simSettings.season + " " ...
        + simSettings.framework + ...
        " T=" + string(T));
    mkdir(path)
end

% Plot simulation behavior for all units
simPlots(path, X, q, sysparams, T, c, lag, printplot);


%% ========================================================================
% SECTION 5: MONTE CARLO SIMS
% ========================================================================

fprintf('Running Monte Carlo Sims.\n');

[V1, V2] = runMonteCarloSims(sysparams, simSettings.bounds, std_hat, X);

fprintf('Simulation complete.\n');

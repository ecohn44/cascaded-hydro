%% Author: Eliza Cohn
% Date: October 2025
% Description: Main driver for syntheitc simulations 
% Paper: Optimization of Cascaded Hydroelectric Systems under DDU

clear; clc; close all;

%% ========================================================================
% SECTION 1: DATA LOADING AND PARAMETER DEFINITION
% ========================================================================

% Static Parameters (PowerProd)
eta = .9;          % efficiency of release-energy conversion
rho_w = 1000;      % density of water [kg/m^3]
g = 9.8;           % acceleration due to gravity [m/s^2]
c = eta*rho_w*g/3.6e9; % power prod coefficient
N = 20;             % number of sub-intervals for piecewise linear approx
n = 1;             % number of cascaded units

% Load inflow data
[inflow, modelparams, sysparams] = dataload(N, n);


%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings
simSettings = initSimSettings("dry", "pwl", "det");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Date range settings 
D = 0;                        % Simulation duration in days
T = 2 + 24*D;                % Number of simulation hours
lag = 1;                      % Number of lag terms in OLS model
year = 2022;                  % Simulation year

% Compute simulation daterange and inflow series
sim_center_date = datetime(year, 1, 1) + days(modelparams.center_day - 1);
start_date = sim_center_date - hours(T/2) - hours(lag);
end_date   = sim_center_date + hours(T/2 - 1);
inflow_s = inflow(inflow.datetime >= start_date & inflow.datetime <= end_date, :);

% Extract historic inflow timeseries [m3/hr]
q_hist = inflow_s.bon_inflow_m3hr; % Input to the first hydrounit

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% ========================================================================

[sol, X, q, std_hat] = optimization(T, N, n, c, q_hist, lag, ...
    simSettings.framework, modelparams, sysparams);

disp(sol.info)

%% ========================================================================
% SECTION 4: PLOTTING
% ========================================================================

% Toggle for creating folder and plotting
make_dir = false;
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

simPlots(path, X, q, q_hist, sysparams, T, c, lag, printplot);


fprintf('Simulation complete.\n');

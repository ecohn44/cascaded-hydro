%% Author: Eliza Cohn
% Date: October 2025
% Description: Main driver for cascaded hydropower simulations 
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

% Load inflow data
[inflow, modelparams, sysparams] = dataload(N);


%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings
simSettings = initSimSettings("dry", "pwl", "ddu");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Date range settings 
D = 1;                        % Simulation duration in days
T = 12 + 24*D;                % Number of simulation hours
lag = 1;                      % Number of lag terms in OLS model
year = 2022;                  % Simulation year

% Compute simulation daterange and inflow series
sim_center_date = datetime(year, 1, 1) + days(modelparams.center_day - 1);
start_date = sim_center_date - hours(T/2) - hours(lag);
end_date   = sim_center_date + hours(T/2 - 1);
inflow_s = inflow(inflow.datetime >= start_date & inflow.datetime <= end_date, :);

% Extract historic inflow timeseries [m3/hr]
q = [inflow_s.bon_inflow_m3hr, inflow_s.tda_inflow_m3hr]; 

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% ========================================================================

[model, obj, X, std_hat] = optimization(T, N, c, q, lag, ...
    simSettings.framework, modelparams, sysparams);


%% ========================================================================
% SECTION 4: PLOTTING
% ========================================================================

% Toggle for creating folder and plotting
make_dir = true;
printplot = true; 

if make_dir
    dir_path = "./plots/";
    stamp = datestr(now,'mm-dd-yyyy HH.MM.SS');
    path = fullfile(dir_path, stamp + " " + simSettings.season + " " ...
        + simSettings.framework + ...
        " T=" + string(T));
    mkdir(path)
end

% Plot simulation behavior for all units
if printplot
    simPlots(path, X, q, sysparams, T, c, lag);
end

fprintf('Simulation complete.\n');

%% Author: Eliza Cohn
% Date: October 2025
% Description: Main driver for cascaded hydropower simulations 
% Paper: Optimization of Cascaded Hydroelectric Systems under DDU

clear; clc; close all;

%% ========================================================================
% SECTION 1: DATA LOADING AND PARAMETER DEFINITION
% ========================================================================

% Load inflow data
[inflow, params, sysparams] = dataload();

% Static Parameters (PowerProd)
eta = .9;          % efficiency of release-energy conversion
rho_w = 1000;      % density of water [kg/m^3]
g = 9.8;           % acceleration due to gravity [m/s^2]


%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

function settings = initSimSettings(season, method, framework)
    validSeasons = ["dry", "wet"];
    validMethods = ["MINLP", "PWL"];
    validFrameworks = ["DET", "DIU", "DDU"];

    if ~ismember(season, validSeasons)
        error('Invalid season. Choose "DRY" or "WET".');
    end
    if ~ismember(method, validMethods)
        error('Invalid method. Choose "MINLP" or "PWL".');
    end
    if ~ismember(framework, validFrameworks)
        error('Invalid framework. Choose "DET", "DIU", or "DDU".');
    end

    settings = struct('season', season, 'method', method, 'framework', framework);
end

% Initialize settings
simSettings = initSimSettings("wet", "PWL", "DET");

% Extract forecasting coefficients 
params = params(strcmp({params.season}, simSettings.season));

% Date range settings 
D = 30;                        % Simulation duration in days
T = 12 + 24*D;                 % Number of simulation hours
lag = 1;                       % Number of lag terms in OLS model
year = 2022;                   % Simulation year

% Compute simulation daterange and inflow series
sim_center_date = datetime(year, 1, 1) + days(params.center_day - 1);
start_date = sim_center_date - hours(T/2) - hours(lag);
end_date   = sim_center_date + hours(T/2 - 1);
inflow_s = inflow(inflow.datetime >= start_date & inflow.datetime <= end_date, :);

% Extract historic inflow timeseries [m3/hr]
q1 = inflow_s.bon_inflow_m3hr;   % Bonneville downstream inflow
q2 = inflow_s.tda_inflow_m3hr;   % Dalles upstream inflow

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% ========================================================================

% model, obj, V1, p1, u1, s1, q1_pred, std_hat, V2, p2, u2, s2 = simulation_loop(q1, q2, framework, method, params)

  %  println("Objective: " * string(obj))

   % println("--- SIMULATION COMPLETE ---")

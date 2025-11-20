%% Author: Eliza Cohn
% Date: November 2025
% Description: Main driver for cascaded hydropower simulations 
% Paper: Optimization of Cascaded Hydroelectric Systems under DDU

tic; 
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
N = 40;             % number of sub-intervals for piecewise linear approx
n = 4;              % number of units in cascaded network 

% Load inflow data 
[modelparams, sysparams, droughtparams] = dataload(n, N);

%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings (season, drought type, linear approximation, uncertainty, bounds)
simSettings = initSimSettings("dry", "extended", "pwl", "ddu", "jcc-bon");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Extract drought simulation mode
droughtparams = droughtparams(strcmp({droughtparams.mode}, simSettings.drought));

% Date range settings 
D = 7;                       % Simulation duration in days
T = 24*D;                     % Number of simulation hours
lag = 3;                      % Number of lag terms in OLS model

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: STREAMFLOW BEHAVIOR
% ========================================================================

% Matrix for local streamflow 
q = zeros(T+lag, n);

% Cascaded drought parameters
baseDrought = droughtparams;
severityScales = [1.0, 0.8, 0.2, 0.1];    
plotDroughtProfiles(baseDrought, severityScales, T, lag, n);

% Simulate drought event 
for i = 1:n
    dp    = baseDrought; 
    scale = severityScales(i);

    if strcmpi(dp.mode, 'extended')
        dp.amp1 = baseDrought.amp1 * scale;
    end

    q(:,i) = droughtSimulator(T, lag, simSettings.season, dp.mode, dp);
end

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
fprintf('Total runtime: %.2f seconds.\n', toc);

function plotDroughtProfiles(baseDrought, severityScales, T, lag, n_units)
% plotDroughtProfiles  Visualize regional extended-drought profile per unit
%
%   baseDrought   : struct like droughtparams(2) for Unit 01
%   severityScales: 1 x n_units vector, e.g. linspace(1.0, 0.5, n_units)
%   T             : main horizon length (in steps, e.g. hours)
%   lag           : extra steps (your inflow lag length)
%   n_units       : number of hydro units

    nTot = T + lag;

    % pull parameters from baseDrought
    daysPerEvent = baseDrought.daysPerEvent;
    tauHours     = baseDrought.tauHours;

    % assume 1 step = 1 hour
    eventLen = max(1, round(daysPerEvent * 24)); % steps per event
    tauSteps = max(1, round(tauHours));         % exponential time constant (steps)

    % time axis
    t = (0:nTot-1)';  % 0-based for convenience

    % Preallocate factors: each column = one unit
    factors = ones(nTot, n_units);

    % For now, assume event starts at t=0 and lasts eventLen steps
    startIdx   = 1;
    fullEndIdx = min(eventLen, nTot);

    tInEvent = (0:(fullEndIdx-startIdx))';  % 0..eventLen-1

    for i = 1:n_units
        amp1_i = baseDrought.amp1 * severityScales(i);

        % exponential shape: 1 -> 1 - amp1_i over time
        decayShape = 1 - exp(- double(tInEvent) / double(tauSteps));
        factors(startIdx:fullEndIdx, i) = 1 - amp1_i * decayShape;
    end

    % ---- Plot ----
    figure('Name', 'Extended Drought Profiles per Unit', 'NumberTitle', 'off');
    tiledlayout(1, n_units, 'Padding', 'compact', 'TileSpacing', 'compact');

    for i = 1:n_units
        nexttile;
        plot(t, factors(:,i), 'LineWidth', 1.5);
        ylim([min(1 - baseDrought.amp1 * severityScales) - 0.05, 1.05]);
        xlim([0, eventLen]);
        xlabel('Time step');
        ylabel('q_{drought} / q_{base}');
        title(sprintf('Unit %02d (amp1=%.2f)', i, baseDrought.amp1 * severityScales(i)));
        grid on;
    end
end
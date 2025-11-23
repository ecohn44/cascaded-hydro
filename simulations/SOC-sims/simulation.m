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
lag = 1;                      % Travel tim ebetween units (hrs)

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: STREAMFLOW BEHAVIOR
% ========================================================================

% Matrix for local streamflow 
q = zeros(T+lag, n);

% Cascaded drought parameters
baseDrought = droughtparams;
severityScales = [1.0, 0.8, 0.2, 0.1];    
plotDroughtProfiles(baseDrought, severityScales, T, lag)

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


function plotDroughtProfiles(baseDrought, severityScales, T, lag)
% plotDroughtProfiles
%   Visualize drought multiplicative factors applied to each unit.
%   This matches the EXTENDED drought logic used in applyExtendedDrought.
%
%   baseDrought     = droughtparams(2)   (extended drought struct)
%   severityScales  = 1 x n_units vector (ex: [1.0 0.7 0.5 0.3])
%

    n_units = numel(severityScales);
    nTot    = T + lag;

    % Extract drought parameters 
    amp1    = baseDrought.amp1;
    nEvents = baseDrought.nEvents;
    daysPerEvent = baseDrought.daysPerEvent;
    tauHours     = baseDrought.tauHours;

    % Convert to steps (1 step = 1 hour)
    eventLen = round(daysPerEvent * 24);
    tauSteps = round(tauHours);
    gapLen   = round((daysPerEvent/2) * 24);   % same gap as simulator

    % Time axis
    t = (1:nTot)';

    % Preallocate drought factors (initially all 1)
    factors = ones(nTot, n_units);

    for e = 1:nEvents

        % event start
        startIdx = (e-1)*(eventLen + gapLen) + 1;
        if startIdx > nTot
            break;
        end

        % event end
        endIdx = min(startIdx + eventLen - 1, nTot);

        % time inside the event: 0 .. (eventLen-1)
        tInEvent = (0:(endIdx - startIdx))';

        % exponential decay: 0 → 1 over the event
        decayShape = 1 - exp(- double(tInEvent) / double(tauSteps));

        % apply drought to each unit
        for u = 1:n_units
            amp_u = amp1 * severityScales(u);
            factors(startIdx:endIdx, u) = 1 - amp_u * decayShape;
        end
    end

    figure('Name','Drought Profiles','NumberTitle','off');
    tiledlayout(1, n_units, 'TileSpacing','compact','Padding','compact');

    ymin = min(factors(:)) - 0.05;

    for u = 1:n_units
        nexttile;
        plot(t, factors(:,u), 'LineWidth', 1.8);
        ylim([ymin 1.05]);
        xlabel('Hour');
        ylabel('q_{drought}/q_{base}');
        title(sprintf('Unit %d (scale=%.2f)', u, severityScales(u)));
        grid on;
    end
end

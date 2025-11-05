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

% Load inflow data 
[modelparams, sysparams] = dataload(N);

%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings (season, linear approximation, uncertainty, bounds)
simSettings = initSimSettings("dry", "pwl", "det", "jcc-bon");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Date range settings 
D = 2;                        % Simulation duration in days
T = 24*D;                     % Number of simulation hours
lag = 1;                      % Number of lag terms in OLS model

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: STREAMFLOW BEHAVIOR
% ========================================================================

q0 = 0.075;    % Dry season base inflow 

% Inflow pulse parameters
amp1 = 0.4;      % 40% drop in inflow
amp2 = 0.3;      % 30% drop in inflow
w1 = 8;          % first drought lasts 8 hours
w2 = 4;          % second drought lasts 4 hours
t0 = [0.3*T, 0.8*T];   % pulses at 30% and 60% of horizon

% q = q0*ones(T+1,1);  % Constant flow 
q = makeInflowPulse(q0, T, lag, t0, amp1, amp2, w1, w2, modelparams.season);

%% ========================================================================
% SECTION 4: OPTIMIZATION FRAMEWORK
% ========================================================================

% [model, obj, X, std_hat, phi_vals, alpha_vals, U_eff] = optimization(T, N, c, q, lag, ...
%     simSettings.framework, simSettings.bounds, modelparams, sysparams);

scale = 1; % Scale safety bounds 
% scale_ddu = 25; % Scale gamma 

[model, obj, X, std_hat, V_eff] = baseOptimization(T, N, c, q, lag, scale, ...
    simSettings.framework, simSettings.bounds, modelparams, sysparams);

% Extract q2 reference inflow
q(:,2) = [0; X(:,3) + X(:,4)];

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
simPlots(path, X, [], q, sysparams, T, c, lag, printplot);

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

%% ========================================================================
% SECTION 6: MONTE CARLO SIMS
% ========================================================================

monte_carlo = true;

if monte_carlo
    fprintf('Running Monte Carlo Sims.\n');
    [V1, V2] = runMonteCarloSims(sysparams, simSettings.bounds, std_hat, X, path, printplot);
end

fprintf('Simulation complete.\n');


function q = makeInflowPulse(q0, T, lag, t0, amp1, amp2, w1, w2, season)
    % makeInflowPulse  Generate inflow with up to two pulses
    % q0     baseline (scalar or vector length T+lag)
    % T,lag  horizon and lag
    % t0     pulse positions (fractions <1 or absolute indices)
    % amp1,amp2  fractional amplitudes (0.2 = +20%)
    % w1,w2      pulse widths (in hours)
    % season     'wet' → positive pulses, 'dry' → negative pulses
    %
    % Returns q (T+lag x 1) >= 0
    
    n = T + lag;
    
    if isscalar(q0)
        q = q0 * ones(n,1);
    else
        q = q0(:);
    end
    
    t0 = t0(:)';  % ensure row vector
    idx = @(x) round((x < 1) .* (lag + x*T) + (x >= 1) .* x);
    t_idx = max(1, min(n, arrayfun(idx, t0)));
    
    sign_mult = 1;
    if strcmpi(season, 'dry')
        sign_mult = -1;
    end
    
    for k = 1:numel(t_idx)
        if k == 1
            amp = sign_mult * max(0, amp1);
            wid = max(0, round(w1));
        else
            amp = sign_mult * max(0, amp2);
            wid = max(0, round(w2));
        end
        if amp == 0 || wid == 0
            continue
        end
        half = floor((wid-1)/2);
        i1 = max(1, t_idx(k) - half);
        i2 = min(n, i1 + wid - 1);
        q(i1:i2) = q(i1:i2) .* (1 + amp);
    end
    
    q = max(q, 0);
    end
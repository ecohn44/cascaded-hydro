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
simSettings = initSimSettings("dry", "extended", "pwl", "ddu", "jcc-ssh");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Extract drought simulation mode
droughtparams = droughtparams(strcmp({droughtparams.mode}, simSettings.drought));

% Date range settings 
D = 3;                       % Simulation duration in days
T = 24*D;                     % Number of simulation hours
lag = 3;                      % Travel time between units (hrs)

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3A: STREAMFLOW BEHAVIOR
% ========================================================================

% Matrix for local streamflow 
q = zeros(T+lag, n);

% Cascaded drought parameters
baseDrought = droughtparams;
severityScales = makeSeverityScales(n);    

% Simulate drought event 
for i = 1:n
    dp    = baseDrought; 
    scale = severityScales(i);

    if strcmpi(dp.mode, 'extended')
        dp.amp1 = baseDrought.amp1 * scale;
    end

    dp.shiftMult = i - 1;

    q(:,i) = droughtSimulator(T, lag, simSettings.season, dp.mode, dp);
end

% Plot streamflow profiles
plotStreamflows(q)

% Offline DIU covariance and std devs
Sigma_diu                 = cov(q);                  % n×n covariance (DIU)
modelparams.Sigma_diu     = Sigma_diu;               % store full covariance
modelparams.sigma_diu     = sqrt(diag(Sigma_diu));   % n×1 std devs

% Offline correlation matrix (for DDU CCC structure)
modelparams.Rcorr         = corr(q);                 % n×n correlation

% Check PD-ness
if isPD(Sigma_diu)
    disp('Σ is PD');
else
    disp('Σ is NOT PD');
end 

%% ========================================================================
% SECTION 4: OPTIMIZATION FRAMEWORK
% ========================================================================

scale = 1; % Scale safety bounds 
% scale_ddu = 25; % Scale gamma 

[model, obj, X, std_hat, V_eff, phi_vals, alpha_vals] = genOptimization(T, N, c, q, lag, scale, ...
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
if simSettings.bounds == "jcc-ssh"
    results_dir = "./resultsSSH/";
else
    results_dir = "./resultsBonferroni/";
end 

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


function plotStreamflows(q)
    % q: T x n matrix, each column is a streamflow time series

    [T, n] = size(q);
    t = (1:T)';              % simple time index; replace with real time if you have it

    figure;
    for i = 1:n
        subplot(n, 1, i);
        plot(t, q(:, i), 'LineWidth', 1.8);
        ylim([0 0.05]); 
        
        ylabel(sprintf('q_%d', i), 'FontSize', 12);
        set(gca, 'FontSize', 12); 

        if i == 1
            title('Streamflow Time Series', 'FontSize', 16);
        end
        if i == n
            xlabel('Time (hour)', 'FontSize', 12);
        else
            set(gca, 'XTickLabel', []);  % hide x labels for middle plots
        end

        grid on;
    end
end


function severityScales = makeSeverityScales(N)

    % Split into two halves
    half = floor(N/2);         % first half length
    rest = N - half;           % second half length
    
    % First half: 1.0, 0.9, 0.8, ..., decreasing by 0.1
    firstHalf  = 1.0 - 0.1*(0:half-1);
    
    % Second half: 0.1, 0.2, 0.3, ..., increasing by 0.1
    secondHalf = 0.1*(1:rest);
    
    % Combine
    severityScales = [firstHalf flip(secondHalf)];
end


function tf = isPD(A)
    % Symmetrize first
    A = (A + A')/2;      
    [~,p] = chol(A);
    tf = (p == 0);        % true if PD
end
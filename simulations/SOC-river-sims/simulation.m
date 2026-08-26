%% Author: Eliza Cohn
% Date: August 2026
% Description: Main driver for cascaded hydropower simulations 

tic; 
clear; clc; close all;

addpath('/Library/gurobi1202/macos_universal2/matlab');
addpath(genpath('/Users/elizacohn/Documents/YALMIP-master'))

% Add shared functions to file path 
thisFilePath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisFilePath, '..', 'functions')));

%% ========================================================================
% SECTION 1: DATA LOADING AND PARAMETER DEFINITION
% ========================================================================

% Toggle for creating folder and plotting
make_dir = true;
printplot = false; 
save_mat = true; 
save_streamflow = false;

% Static parameters 
eta = .9;           % efficiency of release-energy conversion
rho_w = 1000;       % density of water [kg/m^3]
g = 9.8;            % acceleration due to gravity [m/s^2]
c = 1;              % power prod coefficient (c = eta*rho_w*g/3.6e9)
eps = 0.05;         % risk tolerance 

% Load inflow data 
[inflow, soc, modelparams, sysparams] = dataload();

%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings (season, chance constrained solution, uncertainty form)
simSettings = initSimSettings("dry", "det", "det");
 
% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Date range settings            
D = 30;                      % Number of simulation days 
T = D*24;                    % Number of simulation hours
lag = 2;                     % Travel time between units (hrs)
years = [2021,2025];           % Simulation years

% Create path to store results  
if simSettings.bounds == "jcc-ssh"
    results_dir = "./resultsSSH/";
elseif simSettings.bounds == "jcc-bon"
    results_dir = "./resultsBonferroni/";
elseif simSettings.bounds == "det"
    results_dir = "./resultsOracle/";
else 
    warning('Results directory does not exist');
end 

for y = 1:length(years)
    year = years(y);
    fprintf('Running simulation for year: %d\n', year);
    
    stamp = string(year) + "_" + simSettings.season + "_" + simSettings.framework + "_T" + string(T);
    
    % Make plot directory for current simulation run 
    if make_dir
        path = fullfile(results_dir, stamp);
        if ~exist(path, 'dir')
            mkdir(path)
        end
    end
    
    fprintf('Data loading complete.\n');
    
    %% ========================================================================
    % SECTION 3: STREAMFLOW BEHAVIOR
    % ========================================================================
    
    % Compute simulation daterange and input time series
    start_date = datetime(year, 1, 1) + days(modelparams.start_day-1);
    end_date   = start_date + hours(T-1); 
    inflow_s = inflow(inflow.datetime >= start_date & inflow.datetime <= end_date, :);
    soc_s = soc(soc.datetime >= start_date & soc.datetime <= end_date, :);
    
    % Extract and normalize historic inflow timeseries 
    I = [inflow_s.mcn_inflow, inflow_s.jda_inflow, inflow_s.tda_inflow, inflow_s.bon_inflow];
    SOC = [soc_s.mcn_soc, soc_s.jda_soc, soc_s.tda_soc, soc_s.bon_soc];
    
    % Plot streamflow profiles
    plotStreamflows(I);

    % Plot SOC reference trajectories 
    plotSOCs(SOC);
    
    
    %% ========================================================================
    % SECTION 4: OPTIMIZATION FRAMEWORK
    % ========================================================================
    
    thetas = [1.0, 0.1, 0.01];  %Storage tracking penalty
    
    for i = 1:length(thetas)
        theta = thetas(i);
        fprintf('Running simulation for theta: %d\n', theta);

        [model, obj, X] = oracleGurobi(T, c, I', SOC', theta, sysparams);
    
        % Plot simulation behavior for all units
        % simPlots(path, X, SOC, sysparams, T, c, printplot);
        
        % Save results 
        if save_mat
            fname = sprintf('results_theta%d.mat', round(100*theta));
            savepath = results_dir + "/" + stamp + "/";
            save(fullfile(savepath, fname), 'X', 'sysparams', 'T', 'theta');
        end 
        
        
        fprintf('Simulation complete.\n');
        fprintf('Total runtime: %.2f seconds.\n', toc);
    end 
end

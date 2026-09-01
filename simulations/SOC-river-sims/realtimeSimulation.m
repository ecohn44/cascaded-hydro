%% Author: Eliza Cohn
% Date: August 2026
% Description: Main driver for oracle multi-period cascaded hydropower simulations 

tic; 
clear; clc; close all;

addpath('/Library/gurobi1303/macos_universal2/matlab');
addpath(genpath('/Users/elizacohn/Documents/YALMIP-master'))

% Add shared functions to file path 
thisFilePath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisFilePath, '..', 'functions')));

%% ========================================================================
% SECTION 1: DATA LOADING AND PARAMETER DEFINITION
% ========================================================================

% Toggle for creating folder and plotting
make_dir = false;
printplot = false; 
save_mat = false; 
save_streamflow = false;

% Static parameters 
eta = .9;           % efficiency of release-energy conversion
rho_w = 1000;       % density of water [kg/m^3]
g = 9.8;            % acceleration due to gravity [m/s^2]
c = 1;              % power prod coefficient (c = eta*rho_w*g/3.6e9)
eps = 0.05;         % risk tolerance 

% Load inflow data 
[inflow, soc, modelparams, sysparams] = dataload();

% Load reference trajectories 
years = 2018:2024;  
theta = 1;
T = 2160;

results_folder = "resultsOracle";
for k = 1:length(years)
    
    % Define scenario file folder 
    folder = fullfile(results_folder, sprintf('%d_dry_det_T%d',years(k),T));

    % Load in .mat file 
    file = fullfile(folder, sprintf('results_theta%d.mat',round(100*theta)));

    % Unpack simulation results 
    S = load(file,'X','sysparams');

    % Select columns 1, 6, 11, 16 (recorded volumes for each scenario)
    SOC(:,:,k) = S.X(:,1:5:end)';

    % Select columns 5, 10, 15, 20 (historical inflow for each scenario)
    I(:,:,k) = S.X(:,5:5:end)';

end

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
year = 2022;                 % Simulation year

% Create path to store results  
if simSettings.bounds == "jcc-ssh"
    results_dir = "./resultsSSH/";
elseif simSettings.bounds == "jcc-bon"
    results_dir = "./resultsBonferroni/";
elseif simSettings.bounds == "det"
    results_dir = "./resultsRealTime/";
else 
    warning('Results directory does not exist');
end 

fprintf('Running simulation for year: %d\n', year);

% To do: Extract I and SOC for year in year   

% Plot streamflow profiles
plotStreamflows(I);

% Plot SOC reference trajectories 
plotSOCs(SOC);
    
    
%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% ========================================================================

[model, obj, X] = oracleGurobi(T, c, I', SOC', theta, sysparams);

% Plot simulation behavior for all units
% simPlots(path, X, SOC, sysparams, T, c, printplot);

fprintf('Simulation complete.\n');
fprintf('Total runtime: %.2f seconds.\n', toc);

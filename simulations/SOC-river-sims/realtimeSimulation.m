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
theta_ref = 1;  % Select which SOC trajectory to pull
T = 2160;
results_folder = "resultsOracle";
n_units = numel(sysparams);

for k = 1:length(years)
    
    % Define scenario file folder 
    folder = fullfile(results_folder, sprintf('%d_dry_det_T%d',years(k),T));

    % Load in .mat file 
    file = fullfile(folder, sprintf('results_theta%d.mat',round(100*theta_ref)));

    % Unpack simulation results 
    S = load(file,'X','sysparams');

    % Select columns 1, 6, 11, 16 (recorded volumes for each scenario)
    SOC_hist(:,:,k) = S.X(:,1:5:end)';

    % Select columns 5, 10, 15, 20 (historical inflow for each scenario)
    I_hist(:,:,k) = S.X(:,5:5:end)';

end

%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings (season, uncertainty form, chance constrained solution)
simSettings = initSimSettings("dry", "ddu", "jcc-bon");
 
% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Date range settings            
D = 7;                       % Number of simulation days 
T = D*24;                    % Number of simulation hours
lag = 1;                     % Travel time between units (hrs)
year = 2023;                 % Simulation year
theta = 0;                 % Real time tracking coefficient 

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

% Extract mean SOC reference trajectory, inflow, and initial conditions  
year_idx   = find(years == year, 1);
I          = I_hist(:, 1:T, year_idx);
SOC_ref    = mean(SOC_hist(:, 1:T, :),3);
SOC_init   = SOC_hist(:, 1, year_idx);

% Plot streamflow profiles
plotStreamflows(I');

% Plot SOC reference trajectories 
plotSOCs(SOC_ref');
    
    
%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% ========================================================================

% Initialize previous states
V_prev = arrayfun(@(s) s.V0, sysparams(:), 'UniformOutput', true);
u_prev = arrayfun(@(s) s.min_ut, sysparams(:), 'UniformOutput', true);

p_history  = zeros(n_units, T);
u_history  = zeros(n_units, T);
V_history  = zeros(n_units, T);
sp_history = zeros(n_units, T);
q_forecast_history  = zeros(n_units, T);
std_hat    = zeros(n_units, T);

for t = 1:T

    % Calculate upstream release at lagged time step
    if t > lag
        up_release = u_history(:,t-lag) + sp_history(:,t-lag);
    else
        up_release = zeros(n_units,1);
    end

    % Calculate previous forecast error
    if t > 1
        I_prev = I(:,t-1);
        q_error = I_prev - q_forecast_history(:,t-1);
    else
        q_error = zeros(n_units,1);
        I_prev = I(:,1);
    end

    [result, obj, X_t, std_hat(:,t)] = realtimeGurobi(t, c, eps, I_prev, q_error, V_prev, u_prev, ...
        SOC_ref(:,t), theta, lag, up_release, sysparams, modelparams, ...
        simSettings.bounds, simSettings.framework);

    if ismember(result.status, {'OPTIMAL','SUBOPTIMAL','TIME_LIMIT'}) && isfield(result, 'x') && ~isempty(result.x)
        p_history(:,t)  = X_t(:,2);
        u_history(:,t)  = X_t(:,3);
        sp_history(:,t) = X_t(:,4);
        q_forecast_history(:,t)  = X_t(:,5);
        % Volume calculated based on actual inflow not forecasted X_t(:,1)
        % V_history(:,t) = V_prev + I(:,t) - u_history(:,t) - sp_history(:,t); 
        V_history(:,t) = X_t(:,1);
    else
        warning('[t=%d] Solver returned %s. Applying fallback.', t, result.status);
        u_history(:,t)  = u_prev;
        V_history(:,t)  = V_prev + I(:,t) - u_prev;
        for i = 1:n_units
            V_history(i,t) = max(sysparams(i).min_V, min(sysparams(i).max_V, V_history(i,t)));
        end
    end

    V_prev = V_history(:,t);
    u_prev = u_history(:,t);

end

% Store Results 
X = [];
for i = 1:n_units
    X = [X, V_history(i,:)', p_history(i,:)', u_history(i,:)', sp_history(i,:)', q_forecast_history(i,:)'];
end

%% ========================================================================
% SECTION 4: DIAGNOSTICS
% ========================================================================

% Reference tracking error
track_error = zeros(n_units, T);
for i = 1:n_units
    track_error(i,:) = abs(V_history(i,:) - SOC_ref(i,:)) / (sysparams(i).max_V - sysparams(i).min_V);
end

% Volume bounds error
V_min = [sysparams.min_V]';
lower_violation = max(0, V_min - V_history);

total_power = sum(p_history(:));
fprintf('\nSystem Power Generation:          %.2f\n', total_power);
fprintf('Mean normalized tracking error:   %.4f\n',  mean(track_error(:)));


%% ========================================================================
% SECTION 5: PLOTTING AND SAVING
% ========================================================================

simPlots(results_dir, X, SOC_ref, sysparams, T, c, printplot);

if save_mat
    if make_dir && ~exist(results_dir, 'dir'); mkdir(results_dir); end
    save(fullfile(results_dir, sprintf('results_rt_year%d_theta%d.mat', ...
        year, round(100*theta))), 'X', 'sysparams');
    fprintf('Results saved to %s\n', results_dir);
end

fprintf('Simulation complete.\n');
fprintf('Total runtime: %.2f seconds.\n', toc);
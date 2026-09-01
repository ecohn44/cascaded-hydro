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
n_units = numel(sysparams);

for k = 1:length(years)
    
    % Define scenario file folder 
    folder = fullfile(results_folder, sprintf('%d_dry_det_T%d',years(k),T));

    % Load in .mat file 
    file = fullfile(folder, sprintf('results_theta%d.mat',round(100*theta)));

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

% Initialize settings (season, chance constrained solution, uncertainty form)
simSettings = initSimSettings("dry", "det", "det");
 
% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Date range settings            
D = 30;                      % Number of simulation days 
T = 5; %D*24;                    % Number of simulation hours
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

% Extract mean SOC reference trajectory, inflow, and initial conditions  
year_idx   = find(years == year, 1);
I          = I_hist(:, 1:T, year_idx);
SOC_ref    = mean(SOC_hist(:, 1:T, :),3);
SOC_init   = SOC_hist(:, 1, year_idx);

% Plot streamflow profiles
plotStreamflows(I);

% Plot SOC reference trajectories 
plotSOCs(SOC_ref);
    
    
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
q_history  = zeros(n_units, T);

for t = 1:T

    [result, ~, X_t] = realtimeGurobi(t, c, I(:,t), V_prev, u_prev, SOC_ref(:,t), theta, sysparams);

    if ismember(result.status, {'OPTIMAL','SUBOPTIMAL','TIME_LIMIT'}) && isfield(result, 'x') && ~isempty(result.x)
        V_history(:,t)  = X_t(:,1);
        p_history(:,t)  = X_t(:,2);
        u_history(:,t)  = X_t(:,3);
        sp_history(:,t) = X_t(:,4);
        q_history(:,t)  = X_t(:,5);
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
    X = [X, V_history(i,:)', p_history(i,:)', u_history(i,:)', sp_history(i,:)', q_history(i,:)'];
end

%% ========================================================================
% SECTION 4: DIAGNOSTICS
% ========================================================================

p_physical = zeros(n_units, T);
for i = 1:n_units
    p_physical(i,:) = c .* sysparams(i).a .* V_history(i,:).^sysparams(i).b .* u_history(i,:);
end

track_error = zeros(n_units, T);
for i = 1:n_units
    track_error(i,:) = abs(V_history(i,:) - SOC_ref(i,:)) / (sysparams(i).max_V - sysparams(i).min_V);
end

p_error    = p_history - p_physical;
total_phys = sum(p_physical(:));
pct_err    = 100 * (sum(p_history(:)) - total_phys) / max(abs(total_phys), 1e-8);

fprintf('\nSystem Power Generation:          %.2f\n', total_phys);
fprintf('Objective overstatement:          %.2f%%\n', pct_err);
fprintf('Maximum absolute power error:     %.4f\n',  max(abs(p_error(:))));
fprintf('Mean normalised tracking error:   %.4f\n',  mean(track_error(:)));


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
%% Author: Eliza Cohn
% Date: August 2026
% Description: Tuning kernel regression parameters using historical data

clear all; close all;

years = 2018:2024;  
theta = 1;
T = 2160;

results_folder = "resultsOracle";

% Load each year
for k = 1:length(years)
    
    % Define scenario file folder 
    folder = fullfile(results_folder, sprintf('%d_dry_det_T%d',years(k),T));

    % Load in .mat file 
    file = fullfile(folder, sprintf('results_theta%d.mat',round(100*theta)));

    % Unpack simulation results 
    S = load(file,'X','sysparams');

    % Select columns 1, 6, 11, 16 (recorded volumes for each scenario)
    soc_ref(:,:,k) = S.X(:,1:5:end)';

    % Select columns 5, 10, 15, 20 (historical inflow for each scenario)
    inflow_norm(:,:,k) = S.X(:,5:5:end)';

end

% Hyperparamter grid sweep
W_grid = 24; %[24 72]; %[ 168 336];
sigma_grid = 1; %[0.05 0.1 0.25 0.5 1 2 5]; %0.01:0.01:0.10;

% Test kernel regression 
[best_W, best_sigma, mean_rmse, fold_rmse] = kernelRegressionTuning(inflow_norm, soc_ref, W_grid, sigma_grid, years);
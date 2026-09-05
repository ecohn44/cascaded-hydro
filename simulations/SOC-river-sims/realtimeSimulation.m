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

% Static parameters 
eta = .9;           % efficiency of release-energy conversion
rho_w = 1000;       % density of water [kg/m^3]
g = 9.8;            % acceleration due to gravity [m/s^2]
c = 1;              % power prod coefficient (c = eta*rho_w*g/3.6e9)
eps = 0.05;         % risk tolerance 

% Load inflow and system data 
[inflow, soc, modelparams, sysparams] = dataload();
V_min = [sysparams.min_V]';
V_max = [sysparams.max_V]';


% Load reference trajectories 
years = 2018:2025;  
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

% Initialize settings (season, uncertainty form, solution alg, tracking reference)
simSettings = initSimSettings("dry", "ddu", "jcc-bon", "mean");
 
% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Date range settings            
D = 2;                       % Number of simulation days 
T = D*24;                    % Number of simulation hours
lag = 3;                     % Travel time between units (hrs)

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

    
%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% =========================================================================

% Monte Carlo Settings
S = 32;                        % Monte Carlo simulations per year 
kappa = 0.25:0.25:2;           % Forecast error
frameworks = ["diu", "ddu"];   % Uncertainty representation 
thetas = [0, 10, 50, 100];         % Real time tracking coefficient 

% Prepare to save results
M = length(frameworks);
K = length(kappa);
Y = length(years);
H = length(thetas);

% Generate common random draws
rng(1);
R = eye(n_units); % Temp - To replace with AR mean model R estimations 
Z = zeros(n_units,T,S,Y,H);
for y = 1:Y
    for s = 1:S
        Z(:,:,s,y) = mvnrnd(zeros(1,n_units),R,T)';
    end
end

results.V     = nan(n_units,T,Y,H,M,K,S);
results.p     = nan(n_units,T,Y,H,M,K,S);
results.u     = nan(n_units,T,Y,H,M,K,S);
results.sp    = nan(n_units,T,Y,H,M,K,S);
results.q     = nan(n_units,T,Y,H,M,K,S);
results.std   = nan(n_units,T,Y,H,M,K,S);
results.SOC_mean = nan(n_units,T,Y);

% Track infeasible outcomes 
failed = false(Y, H, M, K, S);
failure_time = nan(Y, H, M, K, S);

for y = 1:Y

    year = years(y);
    fprintf('Running simulation for year: %d\n', year);
    
    % Hydrograph and initial conditions
    I          = I_hist(:, 1:T, y);
    SOC_init   = SOC_hist(:, 1, y);

    % Generate reference trajectory from trainin data
    train_years = setdiff(1:Y,y);
    SOC_mean   = mean(SOC_hist(:, 1:T, train_years),3);
    SOC_p10    = prctile(SOC_hist(:, 1:T, train_years),10,3);
    SOC_p90    = prctile(SOC_hist(:, 1:T, train_years),90,3);
    SOC_ref = [SOC_mean; SOC_mean];
    
    for h = 1:H
        theta = thetas(h);
        fprintf('Tracking Coefficient: %d\n', theta);

        for m = 1:M
            framework = frameworks(m);
            fprintf('Uncertainty Framework: %s\n', framework);
            
            for k = 1:K
                fprintf('Forecast Error Level: %d\n', kappa(k));

                for s = 1:S
                    disp(s)
                    % Initialize previous states
                    V_prev = SOC_init; %arrayfun(@(s) s.V0, sysparams(:), 'UniformOutput', true);
                    u_prev  = arrayfun(@(s) s.min_ut, sysparams(:), 'UniformOutput', true);
                    
                    % Initialize histories
                    V_history = nan(n_units,T);
                    p_history = nan(n_units,T);
                    u_history = nan(n_units,T);
                    sp_history = nan(n_units,T);
                    q_mean  = nan(n_units,T);
                    std_hat = nan(n_units,T);
                    q_error = zeros(n_units,T);
                    q_real     = nan(n_units, T);
                                     
                    for t = 1:T
                    
                        % Calculate upstream release at lagged time step
                        if t > lag
                            up_release = u_history(:,t-lag) + sp_history(:,t-lag);
                        else
                            up_release = zeros(n_units,1);
                        end
                    
                        % Calculate previous forecast error
                        if t > 1
                            error_prev = q_error(:,t-1);
                            I_prev = q_real(:,t-1);
                        else
                            error_prev = zeros(n_units,1);
                            I_prev = I(:,1);
                        end
        
                        % Use historical inflow for the first reservoir
                        I_prev(1) = I(1,max(t-1,1));
                                    
                        [result, obj, X_t, std_hat(:,t)] = realtimeGurobi(t, c, eps, I_prev, error_prev, V_prev, u_prev, ...
                            SOC_ref(:,t), theta, lag, up_release, sysparams, modelparams, ...
                            simSettings.bounds, framework, simSettings.ref);
            
                        if ismember(result.status, {'OPTIMAL','SUBOPTIMAL','TIME_LIMIT'}) && isfield(result, 'x') && ~isempty(result.x) 
                            % Calculate realized flow 
                            sigma_ddu = forecast_error(t, error_prev, up_release, "ddu", modelparams, sysparams);
                            q_mean(:,t) = X_t(:,5);
                            q_real(:,t) = max(0,q_mean(:,t) + kappa(k)*sigma_ddu(:).*Z(:,t,s,y));
        
                            p_history(:,t)  = X_t(:,2);
                            u_history(:,t)  = X_t(:,3);
                            sp_history(:,t) = X_t(:,4);
                            V_history(:,t) = V_prev + q_real(:,t) - u_history(:,t) - sp_history(:,t); % volume updated with realized inflow
                        else
                            warning('[t=%d] Solver returned %s. Exiting.', t, result.status);
                            failed(y,h,m,k,s) = true;
                            failure_time(y,h,m,k,s) = t;
                            break
                        end
        
        
                        % Check if still within bounds after real time volume updates 
                        if any(V_history(:,t) < V_min | V_history(:,t) > V_max)
                            warning("Volume out of bounds")
                            failed(y,h,m,k,s) = true;
                            failure_time(y,h,m,k,s) = t;

                            % Water available without going below V_min
                            available = max(V_prev + q_real(:,t) - V_min,0);
                            
                            % Correct turbine release and spill
                            u_history(:,t) = min(u_history(:,t), available);
                            sp_history(:,t) = min(sp_history(:,t), available-u_history(:,t));
                            
                            % Update volume using corrected releases
                            V_history(:,t) = V_prev + q_real(:,t) - u_history(:,t) - sp_history(:,t);

                            % Spill accounting
                            excess = max(V_history(:,t)-V_max,0);
                            sp_history(:,t) = sp_history(:,t) + excess;
                            V_history(:,t) = V_history(:,t) - excess;

                            % Recalculate power using the corrected release
                            for i = 1:n_units
                            
                                V_norm = (V_prev(i)-sysparams(i).min_V) / (sysparams(i).max_V-sysparams(i).min_V);             
                                V_norm = min(1,max(0,V_norm));              
                                head = sysparams(i).min_h + (sysparams(i).max_h-sysparams(i).min_h) * V_norm^sysparams(i).b;
                            
                                p_history(i,t) = c*head*u_history(i,t);
                            end
                        end
                   
                        V_prev = V_history(:,t);
                        u_prev = u_history(:,t);
                        q_error(:,t) = q_real(:,t) - q_mean(:,t);
                    
                    end
        
                    % Store scenario results
                    results.V(:,:,y,h,m,k,s)   = V_history;
                    results.p(:,:,y,h,m,k,s)   = p_history;
                    results.u(:,:,y,h,m,k,s)   = u_history;
                    results.sp(:,:,y,h,m,k,s)  = sp_history;
                    results.q(:,:,y,h,m,k,s)   = q_real;
                    results.std(:,:,y,h,m,k,s) = std_hat;
        
                end
            end
        end
    end
end

% Store final results 
results.failed = failed;
results.failure_time = failure_time;

results.frameworks = frameworks;
results.kappa = kappa;
results.SOC_mean(:,:,y) = SOC_mean;
results.SOC_p10 = SOC_p10;
results.SOC_p90 = SOC_p90;
results.sysparams = sysparams;
results.mean_inflow(y) = mean(I(1,:));

save(fullfile(results_dir,'monteCarloResults.mat'), 'results','-v7.3');

%% ========================================================================
% SECTION 4: DIAGNOSTICS
% ========================================================================

%{

% Plot streamflow profiles
plotStreamflows(I');

% Plot SOC reference trajectories 
plotSOCs(SOC_mean');

% Store Results 
X = [];
for i = 1:n_units
    X = [X, V_history(i,:)', p_history(i,:)', u_history(i,:)', sp_history(i,:)', q_mean(i,:)'];
end

simPlots(results_dir, X, SOC_mean, SOC_p10, SOC_p90, sysparams, T, c, printplot);


total_power = sum(p_history(:));
fprintf('\nSystem Power Generation:          %.2f\n', total_power);
fprintf('Mean normalized tracking error:   %.4f\n',  mean(track_error(:)));


%}

fprintf('Simulation complete.\n');
fprintf('Total runtime: %.2f seconds.\n', toc);
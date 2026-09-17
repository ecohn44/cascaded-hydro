%% Author: Eliza Cohn
% Date: August 2026
% Description: Main driver for real-time cascaded hydropower Monte Carlo
%              simulation using realtimeGurobiMC (McCormick envelope LP).
%              results.p stores PHYSICAL power derived from realized V and u.

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
printplot = false; 
save_mat  = true; 

% Static parameters 
eta   = .9;                     % efficiency of release-energy conversion
rho_w = 1000;                   % density of water [kg/m^3]
g     = 9.8;                    % acceleration due to gravity [m/s^2]
c     =  eta*rho_w*g*28.32/1e6;  % power production coefficient for kcfs
eps   = 0.05;                    % risk tolerance 
n_units = 4;

% Load inflow and system data 
[inflow, soc, modelparams, sysparams] = dataload();
soc_cols = {'mcn_soc','jda_soc','tda_soc','bon_soc'};
V_min = [sysparams.min_V]';
V_max = [sysparams.max_V]';
kV = [sysparams.kV]';

%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings (season, uncertainty form, solution alg, tracking ref)
simSettings = initSimSettings("dry", "ddu", "jcc-ssh", "mean");

% Date range settings            
D   = 2;     % Number of simulation days 
T   = D*24;  % Number of simulation hours
lag = 2;     % Travel time between units (hrs)

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

% Development and final test split
training_years = 2018:2023;
test_year      = 2024;
test_mode      = true;

if test_mode
    years = test_year;
else
    years = training_years;
end

%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% ==========================================                                                                                                                                                                                  ==============================

% Monte Carlo Settings
S          = 5; %10;                    % Monte Carlo simulations per year
kappa      = 1:.5:2;               % Forecast error scaling
frameworks = ["diu","ddu"];       % Uncertainty representation
thetas     = [0, 3, 5]; %0:1:5;    % Real-time tracking coefficient

% Prepare to save results
M = length(frameworks);
K = length(kappa);
Y = length(years);
H = length(thetas);

% Generate common random draws
rng(1);
R = eye(n_units);
Z = zeros(n_units,T,S,Y);
for y = 1:Y
    for s = 1:S
        Z(:,:,s,y) = mvnrnd(zeros(1, n_units), R, T)';
    end
end

% Pre-allocate results arrays
results.V   = nan(n_units, T, Y, H, M, K, S);
results.p   = nan(n_units, T, Y, H, M, K, S);  % physical power
results.u   = nan(n_units, T, Y, H, M, K, S);
results.sp  = nan(n_units, T, Y, H, M, K, S);
results.q   = nan(n_units, T, Y, H, M, K, S);
results.std = nan(n_units, T, Y, H, M, K, S);
results.IVI = nan(n_units,T,Y,H,M,K,S);
results.SOC_mean    = nan(n_units, T, Y);
results.SOC_p10     = nan(n_units, T, Y);
results.SOC_p90     = nan(n_units, T, Y);
results.mean_inflow = nan(1, Y);

% Track infeasible outcomes 
failed       = false(Y, H, M, K, S);
failure_time = nan(Y, H, M, K, S);

for y = 1:Y

    year = years(y);
    fprintf('Running simulation for year: %d\n', year);

    % Get scenario year 
    idx_I   = inflow.data_year == year;
    idx_SOC = soc.data_year == year;

    % Select requested historical year
    I = inflow{idx_I,{'mcn_inflow','jda_inflow','tda_inflow','bon_inflow'}}';
    SOC_year = soc{idx_SOC,{'mcn_soc','jda_soc','tda_soc','bon_soc'}}';
    
    % Time index 
    I = I(:,1:T); 

    % Use all training years
    ref_years = training_years;
    
    % During validation, exclude the year currently being tested
    if ~test_mode
        ref_years = training_years(training_years ~= year);
    end
    
    SOC_train = zeros(4, T, length(ref_years));
    for k = 1:length(ref_years)
        SOC_year = soc{soc.data_year == ref_years(k), soc_cols}';
        SOC_train(:,:,k) = SOC_year(:,1:T);
    end

    SOC_mean = mean(SOC_train, 3);
    SOC_p10  = prctile(SOC_train, 10, 3);
    SOC_p90  = prctile(SOC_train, 90, 3);
    SOC_ref  = [SOC_mean; SOC_mean];

    % Init conditions from tracking 
    SOC_init = SOC_p10(:,1);
    
    results.SOC_mean(:,:,y) = SOC_mean;
    results.SOC_p10(:,:,y)  = SOC_p10;
    results.SOC_p90(:,:,y)  = SOC_p90;
    results.mean_inflow(y)   = mean(I(1,:));

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
                    V_prev = SOC_init;
                    u_prev = arrayfun(@(s) s.min_ut, sysparams(:), 'UniformOutput', true);

                    % Initialize histories
                    V_history  = nan(n_units, T);
                    p_history  = nan(n_units, T);
                    u_history  = nan(n_units, T);
                    sp_history = nan(n_units, T);
                    IVI_history = zeros(n_units,T);
                    q_mean     = nan(n_units, T);
                    std_hat    = nan(n_units, T);
                    q_error    = zeros(n_units, T);
                    q_real     = nan(n_units, T);

                    for t = 1:T

                        % Calculate upstream release at lagged time step
                        if t > lag
                            up_release = u_history(:,t-lag) + sp_history(:,t-lag);
                        else
                            up_release = zeros(n_units, 1);
                        end

                        % Calculate upstream ramp rates at lagged time step 
                        if t > lag + 1
                            up_release_prev = u_history(:,t-lag-1) + sp_history(:,t-lag-1);
                            up_ramp = abs(up_release - up_release_prev);
                        else
                            up_ramp = zeros(n_units,1);
                        end

                        % Calculate previous forecast error
                        if t > 1
                            error_prev = q_error(:,t-1);
                            I_prev     = q_real(:,t-1);
                        else
                            error_prev = zeros(n_units, 1);
                            I_prev     = I(:,1);
                        end

                        % Use historical inflow for the first reservoir
                        I_prev(1) = I(1,t);

                        [result, obj, X_t, std_hat(:,t), phi_val, alpha_vals] = realtimeYALMIP( ...
                            t, c, kappa(k), eps, I_prev, error_prev, V_prev, u_prev, ...
                            SOC_ref(:,t), theta, lag, up_release, up_ramp, sysparams, ...
                            modelparams, simSettings.bounds, framework, simSettings.ref);

                        % Store estimated forecast error under ddu               
                        sigma_ddu   = forecast_error(t, kappa(k), error_prev, up_release, up_ramp, "ddu", modelparams, sysparams);
                        q_mean(:,t) = X_t(:,5);
                        q_real(:,t) = max(0, q_mean(:,t) +  sigma_ddu(:).*Z(:,t,s,y));
      
                        % Use optimized decisions when available
                        if result.problem == 0
                            u_history(:,t)  = X_t(:,3);
                            sp_history(:,t) = X_t(:,4);
                        else
                            warning('Solver failed at t=%d. Holding previous releases.', t);
                        
                            % Initially maintain each unit's previous operating state
                            u_history(:,t)  = u_prev;
                            sp_history(:,t) = zeros(n_units,1);
                        end
                        
                        % Storage produced by those decisions
                        V_raw = V_prev + kV.*(  q_real(:,t) - u_history(:,t) - sp_history(:,t));
                        
                        % Record violations before fallback
                        IVI_history(:,t) = max(V_min - V_raw, 0);
                        
                        % Identify only the units at risk
                        failed_units = V_raw < V_min;
                        
                        % Shut down only those units
                        if any(failed_units)
                            warning('Lower-volume fallback at t=%d for units %s.', ...
                                t, mat2str(find(failed_units)'));
                        
                            u_history(failed_units,t)  = 0;
                            sp_history(failed_units,t) = 0;
                        
                            % Recalculate their storage after shutdown
                            V_raw(failed_units) = V_prev(failed_units) ...
                                + kV(failed_units).*q_real(failed_units,t);
                        end
                        
                        % Correct upper violations using spill
                        extra_spill = max((V_raw - V_max)./kV, 0);
                        
                        sp_history(:,t) = sp_history(:,t) + extra_spill;
                        V_history(:,t)  = V_raw - kV.*extra_spill;
                        
                        % Protect against any remaining numerical deficit
                        V_history(:,t) = max(V_history(:,t), V_min);

                        % Calculate physical power
                        for i = 1:n_units
                            head = sysparams(i).a*V_history(i,t) + sysparams(i).b;
                            p_history(i,t) = c*head*u_history(i,t);
                        end

                        % Advance state
                        V_prev          = V_history(:,t);
                        u_prev          = u_history(:,t);
                        q_error(:,t)    = q_real(:,t) - q_mean(:,t);

                    end  % t loop

                    % Store scenario results
                    results.V(:,:,y,h,m,k,s)   = V_history;
                    results.p(:,:,y,h,m,k,s)   = p_history;
                    results.u(:,:,y,h,m,k,s)   = u_history;
                    results.sp(:,:,y,h,m,k,s)  = sp_history;
                    results.q(:,:,y,h,m,k,s)   = q_real;
                    results.std(:,:,y,h,m,k,s) = std_hat;
                    results.IVI(:,:,y,h,m,k,s) = IVI_history;

                    if printplot
                        X = [];
                        for i = 1:n_units
                            X = [X, V_history(i,:)', p_history(i,:)', u_history(i,:)', sp_history(i,:)', q_real(i,:)'];
                        end
                        simPlots(results_dir, X, SOC_mean, SOC_p10, SOC_p90, sysparams, T, c, std_hat, eps, false);
                    end 

                end  % s loop
            end  % k loop
        end  % m loop
    end  % h loop
end  % y loop

% Store final metadata
results.failed       = failed;
results.failure_time = failure_time;
results.frameworks   = frameworks;
results.kappa        = kappa;
results.sysparams    = sysparams;
results.thetas       = thetas;
results.years        = years;

save(fullfile(results_dir, 'monteCarloResults.mat'), 'results', '-v7.3');

%% ========================================================================
% SECTION 4: DIAGNOSTICS
% ========================================================================
%{
plotStreamflows(I');
plotSOCs(SOC_mean');

X = [];
for i = 1:n_units
    X = [X, V_history(i,:)', p_history(i,:)', u_history(i,:)', sp_history(i,:)', q_mean(i,:)'];
end
simPlots(results_dir, X, SOC_mean, SOC_p10, SOC_p90, sysparams, T, c, std_hat, eps, printplot);

total_power = sum(p_history(:));
fprintf('\nSystem Power Generation:          %.2f\n', total_power);
fprintf('Mean normalized tracking error:   %.4f\n',  mean(track_error(:)));

soc_ref = [SOC_p10, SOC_p90]

%}

fprintf('Simulation complete.\n');
fprintf('Total runtime: %.2f seconds.\n', toc);
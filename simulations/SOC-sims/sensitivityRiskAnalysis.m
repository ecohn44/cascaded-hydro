%% ========================================================================
% Author: Eliza Cohn
% Date: March 2026
% Description: Probability Parameter Generator for MC Sims 
% =========================================================================

tic; 
clear; clc; close all;
% rng(0, 'twister');

addpath('/Library/gurobi1202/macos_universal2/matlab');
addpath(genpath('/Users/elizacohn/Documents/YALMIP-master'))

% Add shared functions to file path 
thisFilePath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisFilePath, '..', 'functions')));

%% STEP 1: DATA LOADING AND PARAMETER DEFINITION

% Static parameters 
eta = .9;           % efficiency of release-energy conversion
rho_w = 1000;       % density of water [kg/m^3]
g = 9.8;            % acceleration due to gravity [m/s^2]
c = 1;              % power prod coefficient (c = eta*rho_w*g/3.6e9)
N = 40;             % number of sub-intervals for piecewise linear approx
n = 3;              % number of units in cascaded network 
eps = 0.05;         % risk tolerance 
K = 100;            % number of MC sims

% Load simulation parameters
[modelparams, sysparams, seasonparams] = dataload(n, N);

% Initialize settings (season, drought type, lin a, pprox, uncertainty, sln alg, volume price)
simSettings = initSimSettings("dry", "extended", "pwl", "ddu", "jcc-ssh", "none");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Extract event scenario simulation mode
seasonparams = seasonparams(strcmp({seasonparams.mode}, simSettings.scenario));

% Set baseline flow based on season
seasonparams.q0 = 0.85*sysparams(1).max_ut;
[sysparams.V0] = deal(0.25);

% Date range settings            
D = 3.5;        % Simulation duration in days
T = 24*D;       % Number of simulation hours
lag = 1;        % Travel time between units (hrs)

% Load price data
LMP = ones(T, n); % simulatePrice(T, n, true);

% Struct to hold outputs for each scenario
alpha_set = .1:0.05:0.6;
duration_set = 0.8:.1:2;
phi_sweep = zeros(length(alpha_set), length(duration_set));

fprintf('Data loading complete.\n');

%% Step 2: SAMPLE STREAMFLOW BEHAVIOR & OPTIMIZE
for a = 1:length(alpha_set)
    alpha = alpha_set(a);
    
    for d = 1:length(duration_set)
        eventDur = duration_set(d);

        fprintf('Testing values of alpha = %.2f, D = %.2f\n', alpha, eventDur);
        
        % Matrix for local streamflow 
        q = zeros(T+lag, n);
    
        % Build drought event 
        seasonparams.amp1 = alpha;
        seasonparams.daysPerEvent = eventDur;
        
        % Cascaded drought parameters
        baseStreamflow = seasonparams;
        
        % Simulate cascaded drought events 
        for i = 1:n
            dp    = baseStreamflow; 
            dp.shiftMult = i - 1;
            q(:,i) = scenarioSimulator(T, lag, simSettings.season, dp.mode, dp);
        end
    
        % Offline DIU covariance and std devs
        sigma_diu                 = modelparams.AR_std*ones(n,1);   % n×1 std devs
        modelparams.sigma_diu     = sigma_diu;        
        modelparams.Sigma_diu     = diag(sigma_diu.^2);             % n×n covariance (DIU)
        
        % Offline correlation matrix
        modelparams.Rcorr         = estimateR(T, n, lag, q, modelparams); % n×n correlation
    
        [model, obj, X, std_hat, V_eff, phi_vals, alpha_vals] = genOptimization(T, N, c, q, LMP, lag, 1, ...
        simSettings.framework, simSettings.bounds, modelparams, sysparams, eps, simSettings.volPrice);

        % Manually adjust phi(t = 1) since under init conditions 
        phi_vals(1) = 1;
    
        % Store risk attribution values
        phi_sweep(a,d) =  min(phi_vals);

    end 

end 

% plotRiskSpread(alphas_MC)

% Plot Heat MaP
figure;
imagesc(alpha_set, duration_set, phi_sweep)
set(gca,'YDir','normal');

cblabel = 'Minimum Reliability Guarantee';
xlabel('Drought Intensity (\alpha)', 'FontSize', 12);
ylabel('Drought Duration (D)', 'FontSize', 12);

cb = colorbar;
ylabel(cb, cblabel, 'FontSize', 14);

grid on;
set(gca,'GridColor','k','GridAlpha',0.4,'LineWidth',1.0);


save('phi_sweep.mat', 'phi_sweep');
fprintf('Total runtime: %.2f seconds.\n', toc);

function R = estimateR(T, n, lag, q, m)

    alpha0 = m.AR_const;
    alpha1 = m.AR_coef;

    % Preallocate forecast and residuals
    q_hat = nan(T, n);
    E     = nan(T, n);
    
    % AR(1) forecast (q already starts at lag)
    for t = 1:T
        % Forecast streamflow based on lagged observation 
        q_hat(t,:) = alpha0 + alpha1*q(t,:);

        % Calculate forecast error based on actual observation 
        E(t,:)     = q(t+lag,:) - q_hat(t,:);
    end

    % Compute correlation of forecast errors
    R = corr(E, 'Rows', 'complete');

end
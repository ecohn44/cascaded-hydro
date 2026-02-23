%% Author: Eliza Cohn
% Date: February 2026
% Description: Sensitivity analysis for cascaded hydropower simulations 
% Paper: EXPLICIT RISK ALLOCATION FOR CASCADED HYDROELECTRIC SYSTEMS UNDER EXTREME EVENTS

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

% Static parameters 
eta = .9;           % efficiency of release-energy conversion
rho_w = 1000;       % density of water [kg/m^3]
g = 9.8;            % acceleration due to gravity [m/s^2]
c = 1;              % power prod coefficient (c = eta*rho_w*g/3.6e9)
N = 40;             % number of sub-intervals for piecewise linear approx
n = 3;              % number of units in cascaded network 
eps = 0.05;         % risk tolerance 

% Load inflow data 
[modelparams, sysparams, seasonparams] = dataload(n, N);

%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings (season, drought type, lin approx, uncertainty, sln alg, volume price)
simSettings = initSimSettings("wet", "pulse", "pwl", "ddu", "jcc-bon", "none");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Extract event scenario simulation mode
seasonparams = seasonparams(strcmp({seasonparams.mode}, simSettings.scenario));

% Set baseline flow based on season
if simSettings.season == "dry"
    seasonparams.q0 = 0.9*sysparams(1).max_ut;
    [sysparams.V0] = deal(0.1);
elseif simSettings.season == "wet"
    seasonparams.q0 = 0.95*sysparams(1).max_ut;
    [sysparams.V0] = deal(0.8);
else
    seasonparams.q0 = sysparams(1).max_ut;
end

% Date range settings 
D = 3.5;                      % Simulation duration in days
T = 24*D;                     % Number of simulation hours
lag = 1;                      % Travel time between units (hrs)

% Load price data
LMP = ones(T, n); 

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3A: STREAMFLOW BEHAVIOR
% ========================================================================

% Matrix for local streamflow 
q = zeros(T+lag, n);

% Cascaded drought parameters
baseStreamflow = seasonparams;

% Simulate scaled events 
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

% Check PD-ness for covariance matrix 
if isPD(modelparams.Sigma_diu)
    disp('Σ is PD');
else
    warning('Σ is NOT PD');
end 

% Check PD-ness for correlation matrix 
if isPD(modelparams.Rcorr)
    disp('R is PD');
else
    warning('R is NOT PD');
end 

%% ========================================================================
% SECTION 4: SENSITIVITY ANALYSIS & OPTIMIZATION FRAMEWORK
% ========================================================================

% Record initial parameters
gamma0 = modelparams.gamma;
xi0 = modelparams.alpha; 

% Multiplicative sweep factors 
xi_factors    = [0, 1, 5, 10, 50];
gamma_factors = [0, 0.5, 1, 1.5, 2, 3];

% Build sweep grids
xi_vals    = xi0    .* xi_factors;
gamma_vals = gamma0 .* gamma_factors;
n_xi    = length(xi_vals);
n_gamma = length(gamma_vals);

% Intialize array to store objective values 
J_vals = zeros(n_xi, n_gamma); 
IVI_vals = zeros(n_xi, n_gamma); 

for i = 1:n_xi
    xi = xi_vals(i);
    
    for j = 1:n_gamma
        gamma = gamma_vals(j);

        fprintf('Running xi = %.4f, gamma = %.4f\n', xi, gamma);
        
        % Update model parameters
        modelparams.alpha = xi;
        modelparams.gamma = gamma;
        
        % Run optimization
        [model, obj, X, std_hat, V_eff, phi_vals, alpha_vals] = genOptimization(T, N, c, q, LMP, lag, 1, ...
            simSettings.framework, simSettings.bounds, modelparams, sysparams, eps, simSettings.volPrice);

        % Store objective value
        J_vals(i,j) = obj;

        % Load in historical streamflow 
        for k = 1:n
            base = 5*(k-1);
            X(:, base+5) = q(1+lag:T+lag,k);
        end

        % Run Policy Test Sims 
        [~, ~, ~, IVI] = runPolicyTestSims(sysparams, simSettings.bounds, X, 'DDU', std_hat);

        % Store system-wide constraint violation 
        IVI_vals(i,j) = sum(IVI);
    end
end


%% ========================================================================
% SECTION 5: PLOTTING
% ========================================================================

plotHeat(gamma_vals, xi_vals, J_vals, gamma0, xi0, "obj")
plotHeat(gamma_vals, xi_vals, IVI_vals, gamma0, xi0, "IVI")

fprintf('Simulation complete.\n');
fprintf('Total runtime: %.2f seconds.\n', toc);

function tf = isPD(A)
    % Symmetrize first
    A = (A + A')/2;      
    [~,p] = chol(A);
    tf = (p == 0);        % true if PD
end

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

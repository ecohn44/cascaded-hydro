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
simSettings = initSimSettings("dry", "extended", "pwl", "ddu", "jcc-ssh", "none");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));

% Extract event scenario simulation mode
seasonparams = seasonparams(strcmp({seasonparams.mode}, simSettings.scenario));

% Set baseline flow based on season
if simSettings.season == "dry"
    seasonparams.q0 = 0.85*sysparams(1).max_ut;
    [sysparams.V0] = deal(0.25);
    D = 3.5;  

    seasonparams.amp1 = .35;
end

% Date range settings 
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
omega0 = modelparams.omega; 

% Multiplicative sweep factors 
xi_factors    = linspace(1, 100, 5); 
gamma_factors = linspace(1, 10, 5); 
omega_factors = [0, 1, 2, 5];

% Build sweep grids
xi_vals    = xi0    .* xi_factors;
gamma_vals = gamma0 .* gamma_factors;
omega_vals = omega0 .* omega_factors; 
n_xi    = length(xi_vals);
n_gamma = length(gamma_vals);
n_omega = length(omega_vals);

% Intialize array to store objective values 
J_vals = zeros(n_xi, n_gamma); 
IVI_vals = zeros(n_xi, n_gamma); 

outDir = fullfile(thisFilePath, 'mats');

S1 = load(fullfile(outDir, 'DDU_J_vals.mat'), 'J_vals');
S2 = load(fullfile(outDir, 'DDU_IVI_vals.mat'), 'IVI_vals');
J_vals = S1.J_vals;
IVI_vals = S2.IVI_vals;

%{
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
        [~, ~, ~, IVI] = runPolicyTestSims(sysparams, simSettings.bounds, X, 'DDU');

        % Store system-wide constraint violation 
        IVI_vals(i,j) = sum(IVI);
    end
end
%}

%% ========================================================================
% SECTION 5: PLOTTING
% ========================================================================

xlab = 'Upstream Release Coefficient (\gamma)';
ylab = 'Previous Forecast Error (\xi)';

%plotHeat(gamma_vals, xi_vals, J_vals, gamma0, xi0, "obj", xlab, ylab)
%plotHeat(gamma_vals, xi_vals, rescale(IVI_vals, "volume"), gamma0, xi0, "IVI", xlab, ylab)

% Save results
save(fullfile(outDir, 'DDU_J_vals.mat'), 'J_vals');
save(fullfile(outDir, 'DDU_IVI_vals.mat'), 'IVI_vals');

fprintf('Simulation complete.\n');
fprintf('Total runtime: %.2f seconds.\n', toc);
    

% Average J_vals and IVI_vals over xi (rows)
J_avg = mean(J_vals, 1);      % 1 × n_gamma
IVI_avg = mean(rescale(IVI_vals, "volume"), 1);  % 1 × n_gamma

figure('Color','w'); % white background

markersize = 8;

%% Plot: Mean Objective (J)
subplot(2,1,1); hold on; grid on;
plot(gamma_vals, J_avg, '-', 'Color', [0.85, 0.33, 0.1], 'LineWidth',2, ...   % Dark green
     'Marker','o', 'MarkerFaceColor', [0.85, 0.33, 0.1], 'MarkerSize',markersize);
xlabel(xlab, 'FontSize',14);
ylabel('Mean Generation Increase (%)', 'FontSize',14);
set(gca, 'FontSize',12, 'LineWidth',1.2, 'Box','on');
ylim([15 23])
xlim([min(gamma_vals) max(gamma_vals)]);

%% Plot: Mean IVI
subplot(2,1,2); hold on; grid on;
plot(gamma_vals, IVI_avg/1e5, '-', 'Color', [0.4, 0.4, 0.4], 'LineWidth',2, ...  % Dark blue
     'Marker','o', 'MarkerFaceColor', [0.4, 0.4, 0.4], 'MarkerSize',markersize);
xlabel(xlab, 'FontSize',14);
ylabel('Mean IVI (10^5 m^3)', 'FontSize',14);
set(gca, 'FontSize',12, 'LineWidth',1.2, 'Box','on');
xlim([min(gamma_vals) max(gamma_vals)]);
ylim([5 35])
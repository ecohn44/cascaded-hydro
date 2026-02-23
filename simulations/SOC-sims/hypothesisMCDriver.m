%% Monte Carlo Driver and Benchmarking
% Author: Eliza Cohn
% Description:
%   Benchmark DET (M1), DIU (M2), and DDU (M3) policies for a 4-unit
%   cascade under their corresponding uncertainty frameworks using
%   Monte Carlo simulations. Then report total generation and violation
%   statistics for your table.

clc; clear; close all;

% Add shared functions to file path 
thisFilePath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisFilePath, '..', 'functions')));

%% Simulation settings

mode = 1; % Toggle between testing hypothesis 1 or 2

season = "wet";
baseFolder = './resultsSSH/' + season;
simSettings.bounds = "jcc-bon";
printplot = true;
hist_flood = true; 
path      = "";

n = 3;
N = 40;

% Load in historical streamflow data 
q_flood =load("floodFlow.mat").q_save;

% Load model + system parameters 
[params, sysparams, ~] = dataload(n, N);

% Set baseline flow based on season
if season == "dry"
    [sysparams.V0] = deal(0.1);
elseif season == "wet"
    [sysparams.V0] = deal(0.8);
end

% Policy codes (used in filenames) and labels (M1/M2/M3 for reporting)
polCodes = {'det','diu','ddu'};   % DET, DIU, DDU
polNames = {'M1','M2','M3'};      % labels for results
n_pols = length(polNames);

% Struct to hold MC outputs for each policy
MC = struct();

% For H2 set all variance to DDU estimation 
if mode == 2
    polCode = 'ddu';
    fname = fullfile(baseFolder, sprintf('results_unit1_%s.mat', polCode));
    S = load(fname);          % expects: S.X (T x 5n), S.std_hat (T x n)
    std_all = S.std_hat; 
end

%% Loop over policies and run MC sims
for k = 1:numel(polCodes)
    polCode = polCodes{k};    % 'det', 'diu', 'ddu'
    polName = polNames{k};    % 'M1', 'M2', 'M3'
    fprintf('%s\n', polName);

    policyLabel = sprintf('%s (%s)', polName, upper(polCode)); 

    % Load policies derived under optimization 
    fname = fullfile(baseFolder, sprintf('results_unit1_%s.mat', polCode));
    S = load(fname);          % expects: S.X (T x 5n), S.std_hat (T x n)

    T       = S.T;
    X_all   = S.X;   
    tt      = (1:T)';

    plotInflowMismatchUnit1(X_all, q_flood, policyLabel);

    % load in base streamflow data (no predictions)
    if hist_flood
        for i = 1:n
            base = 5*(i-1);
            X_all(:, base+5) = q_flood(:,i);
        end

    end

    % Load corresponding std for H1
    if mode == 1
        std_all = S.std_hat; 
    end 

    % Run Policy Test Sims 
    [V_sim, u_sim, p_sim, IVI] = runPolicyTestSims(sysparams, simSettings.bounds, X_all, policyLabel, std_all);

    V_plan = X_all(:,1:5:end);
    dV = V_plan - V_sim;

    % Run MC Test Sims
    % [V_sim, u_sim, p_sim, MFV, RLR, IVI, spillStats] = runMonteCarloSims(sysparams, simSettings.bounds, std_all, X_all, path, printplot, policyLabel);

    % Mean power per unit over MC runs 
    p_mean = mean(p_sim, 3); % T x n_units

    % Store in MC struct
    MC.(polName).T       = T;
    MC.(polName).tt      = tt;
    MC.(polName).X       = X_all;
    MC.(polName).std_hat = std_all;
    MC.(polName).V       = V_sim;       % T x n_units x nSim
    MC.(polName).U       = u_sim;       % T x n_units x nSim
    MC.(polName).P       = p_sim;       % T x n_units x nSim
    MC.(polName).P_mean  = p_mean;      % T x n_units
    MC.(polName).MFV     = MFV;         % 1 x n_units 
    MC.(polName).RLR     = RLR;         % 1 x n_units
    MC.(polName).IVI     = IVI;         % 1 x n_units
    fprintf('%s complete.\n\n', polName);
end

%% Benchmarking / Reporting
% Compute:
%   1. Total generation per policy (sum over all units and time)
%   2. Generation per unit
%   3. Mean MFV over all units 
%   4. Mean RLR over all units
%   5. Mean IVI over all units 

totalGen  = zeros(1, n_pols);
unitGen   = zeros(n, n_pols);
meanMFV   = zeros(1, n_pols);
meanRLR   = zeros(1, n_pols);
meanIVI   = zeros(1, n_pols);

for k = 1:n_pols
    
    name = polNames{k};
    P      = MC.(name).P_mean;   % T x n_units
    MFV    = MC.(name).MFV;      % 1 x n_units
    RLR    = MC.(name).RLR;      % 1 x n_units
    IVI    = MC.(name).IVI;      % 1 x n_units

    % Total generation per policy
    totalGen(k)  = sum(P(:));
    unitGen(:,k) = sum(P, 1).';

    % Average metric per policy
    meanMFV(k)  = mean(MFV);
    meanRLR(k)  = mean(RLR);
    meanIVI(k)  = mean(IVI);
end


%% Print summary for table
fprintf('\nMonte Carlo Benchmark Summary\n\n');

fprintf('Mean Integrated Violation Index (IVI):\n');
for k = 1:n_pols
    fprintf('  %s: %.4f\n', polNames{k}, 100*meanIVI(k));
end
fprintf('\n');



function plotInflowMismatchUnit1(X, q_hist, policyLabel)

    % X:        T x (5n) optimization output
    % q_hist:   T x n historical inflow used in replay
    % policyLabel: string for figure title

    T  = size(X,1);
    tt = (1:T)';

    % Extract predicted inflow for Unit 1
    q_pred = X(:,5);        % column 5 = q_1(t)
    q_hist1 = q_hist(:,1);

    dq = q_pred - q_hist1;

    figure('Name','Inflow Comparison Unit 1','NumberTitle','off');
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    % --- Top: Historical vs Predicted ---
    nexttile; hold on; grid on;
    plot(tt, q_hist1, 'k-',  'LineWidth', 2,   'DisplayName','Historical');
    plot(tt, q_pred,  'b--', 'LineWidth', 1.8, 'DisplayName','Predicted');

    ylabel('q_1');
    title('Unit 1 Inflow: Historical vs Predicted');
    legend('Location','best');

    % --- Bottom: Difference ---
    nexttile; hold on; grid on;
    plot(tt, dq, 'r-', 'LineWidth', 2);
    yline(0,'k--','HandleVisibility','off');

    xlabel('Time');
    ylabel('\Delta q_1');
    title('Prediction Error: q_{pred} - q_{hist}');

    sgtitle(sprintf('Forecast vs Historical Inflow — %s', policyLabel));

end

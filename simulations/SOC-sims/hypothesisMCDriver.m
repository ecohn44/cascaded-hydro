%% Hypothesis 01 Monte Carlo Driver and Benchmarking
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

mode = 2; % Toggle between testing hypothesis 1 or 2

baseFolder = './resultsSSH';
simSettings.bounds = "jcc-ssh";
printplot = true;
path      = "";

n = 3;
N = 40;

% Load model + system parameters 
[params, sysparams, ~] = dataload(n, N);

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
    % Load corresponding std for H1
    if mode == 1
        std_all = S.std_hat; 
    end 

    % Run Monte Carlo sims (general n-unit version)
    [V_sim, u_sim, p_sim, MFV, RLR, IVI] = ...
        runMonteCarloSims(sysparams, simSettings.bounds, std_all, X_all, path, printplot, policyLabel);

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
    MC.(polName).IVI     = IVI;  
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

fprintf('Mean Frequency of Violations (MFV):\n');
for k = 1:n_pols
    fprintf('  %s: %.2f%% of runs\n', ...
            polNames{k}, 100*meanMFV(k));
end
fprintf('\n');

fprintf('Mean Run-Level Risk (RLR):\n');
for k = 1:n_pols
    fprintf('  %s: %.4f\n', polNames{k}, 100*meanRLR(k));
end
fprintf('\n');

fprintf('Mean Integrated Violation Index (IVI):\n');
for k = 1:n_pols
    fprintf('  %s: %.4f\n', polNames{k}, 100*meanIVI(k));
end
fprintf('\n');


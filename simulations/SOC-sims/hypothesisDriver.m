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

mode = 1; % Toggle between testing hypothesis 1 or 2

simSettings.bounds = "jcc-bon";
printplot = true;
path      = "";

n = 4;
N = 40;

% Load model + system parameters 
[params, sysparams, ~] = dataload(n, N);

% Policy codes (used in filenames) and labels (M1/M2/M3 for reporting)
polCodes = {'det','diu','ddu'};   % DET, DIU, DDU
polNames = {'M1','M2','M3'};      % labels for results
n_pols = length(polNames);

baseFolder = './results';

% Struct to hold MC outputs for each policy
MC = struct();

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
    std_all = S.std_hat;      
    tt      = (1:T)';

    % Run Monte Carlo sims (general n-unit version)
    [V_sim, u_sim, p_sim, viol_frac, prob_run_violation] = ...
        runMonteCarloSims(sysparams, simSettings.bounds, std_all, X_all, path, printplot, policyLabel);

    % Mean power per unit over MC runs (T x n)
    p_mean = mean(p_sim, 3);

    % Store in MC struct
    MC.(polName).T       = T;
    MC.(polName).tt      = tt;
    MC.(polName).X       = X_all;
    MC.(polName).std_hat = std_all;
    MC.(polName).V       = V_sim;       % T x n_units x nSim
    MC.(polName).U       = u_sim;       % T x n_units x nSim
    MC.(polName).P       = p_sim;       % T x n_units x nSim
    MC.(polName).P_mean  = p_mean;      % T x n_units
    MC.(polName).D_time  = viol_frac;   % T x n_units (fraction of runs violating)
    MC.(polName).D_run  = prob_run_violation; % 1 x n (run-level)

    fprintf('%s complete.\n\n', polName);
end

%% Benchmarking / Reporting
% Compute:
%   - Total generation per policy (sum over all units and time)
%   - Generation per unit
%   - Total violation counts per policy (sum of fractions over time)
%   - Violation counts per unit

totalGen  = zeros(1, n_pols);
unitGen   = zeros(n, n_pols);
totalViol = zeros(1, n_pols);
unitViol  = zeros(n, n_pols);
runViolProb = zeros(n, n_pols);

for k = 1:n_pols
    name = polNames{k};

    P      = MC.(name).P_mean;   % T x n_units
    D_time = MC.(name).D_time;   % T x n_units
    D_run  = MC.(name).D_run;    % 1 x n_units

    % Total generation
    totalGen(k)  = sum(P(:));
    unitGen(:,k) = sum(P, 1).';

    totalViol_time(k)  = sum(D_time(:));
    unitViol_time(:,k) = sum(D_time, 1).';
    runViolProb(:,k) = D_run.';   % n_units x n_pols
end


%% Print summary for table
fprintf('\nMonte Carlo Benchmark Summary (Hypothesis 01)\n\n');
fprintf('Total Generation (p.u.):\n');
for k = 1:n_pols
    fprintf('  %s : %.4f\n', polNames{k}, totalGen(k));
end
fprintf('\n');

fprintf('Generation per Unit (p.u.):\n');
for i = 1:n
    fprintf('  Unit %d:\n', i);
    for k = 1:n_pols
        fprintf('    %s : %.4f\n', polNames{k}, unitGen(i,k));
    end
end
fprintf('\n');

fprintf('Run-Level Violation Probability:\n');
fprintf('(Probability a run violates at least once)\n');
for i = 1:n
    fprintf('  Unit %d:\n', i);
    for k = 1:n_pols
        fprintf('    %s : %.2f%%%% of runs\n', ...
                polNames{k}, 100*runViolProb(i,k));
    end
end
fprintf('\n');

fprintf('Total Violation (Metric B: sum of per-time-step fractions):\n');
for k = 1:n_pols
    fprintf('  %s : %.4f\n', polNames{k}, totalViol_time(k));
end
fprintf('\n');

fprintf('Violation (Metric B) per Unit (sum of fractions over time):\n');
for i = 1:n
    fprintf('  Unit %d:\n', i);
    for k = 1:n_pols
        fprintf('    %s : %.4f\n', polNames{k}, unitViol_time(i,k));
    end
end
fprintf('\n');

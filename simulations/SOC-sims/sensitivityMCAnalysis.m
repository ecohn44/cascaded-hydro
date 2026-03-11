%% ========================================================================
% Author: Eliza Cohn
% Date: March 2026
% Description: Probability Parameter Generator for MC Sims 
% =========================================================================

tic; 
clear; clc; close all;
rng(0, 'twister');

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
K = 20000;            % number of MC sims

% Load simulation parameters
[modelparams, sysparams, seasonparams] = dataload(n, N);

% Policy test settings
season = "dry";
scenario = "extended";
bounds = "jcc-bon";
baseFolder = './resultsBonferroni/' + season;

% Policy codes (used in filenames) 
polCodes = {'det','diu','ddu'};   
polNames = {'M1','M2','M3'};     
n_pols = length(polNames);

% Define a struct for policy loading 
policies = struct();

% Load optimal policy 
for m = 1:numel(polCodes)
    
    polCode = polCodes{m};    
    polName = polNames{m};    
    policyLabel = sprintf('%s (%s)', polName, upper(polCode)); 

    % Load policies derived under optimization 
    fname = fullfile(baseFolder, sprintf('results_unit1_%s.mat', polCode));
    S = load(fname);   

    % Store optimal policies
    policies(m).code  = polCode;
    policies(m).name  = polName;
    policies(m).X     = S.X;
end 

% Struct to hold MC outputs for each policy
IVI_MC = zeros(K, n_pols);
P_MC = zeros(K, n_pols);
viol_MC = zeros(K, n_pols);

q0_MC = zeros(K, 1);
amp_MC = zeros(K, 1);
D_MC = zeros(K, 1);

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, season));

% Extract event scenario simulation mode
seasonparams = seasonparams(strcmp({seasonparams.mode}, scenario));

% Set SoC and horizon length 
[sysparams.V0] = deal(0.25);
D = 3.5;

% Date range settings 
T = 24*D;                     % Number of simulation hours
lag = 1;                      % Travel time between units (hrs)

% Load price data
LMP = ones(T, n); 

fprintf('Data loading complete.\n');

%% Step 2: SAMPLE STREAMFLOW BEHAVIOR
for k = 1:K

    % Matrix for local streamflow 
    q = zeros(T+lag, n);

    seasonparams.q0 = 0.0425 + 0.004*randn();
    seasonparams.amp1 = betarnd(1.5, 8.5);
    seasonparams.daysPerEvent = gamrnd(4, 0.25);

    q0_MC(k) = seasonparams.q0;
    amp_MC(k) = seasonparams.amp1;
    D_MC(k) = seasonparams.daysPerEvent;
    
    % Cascaded drought parameters
    baseStreamflow = seasonparams;
    
    % Simulate cascaded drought events 
    for i = 1:n
        dp    = baseStreamflow; 
        dp.shiftMult = i - 1;
        q(:,i) = scenarioSimulator(T, lag, season, dp.mode, dp);
    end
    
    
    %% Step 3: TEST OPTIMAL POLICIES
    % Loop over Policies and run MC sims
    for m = 1:n_pols
       
        % Load optimal policy 
        X_opt   = policies(m).X ;   
    
        % Load in sampled streamflow data 
        for i = 1:n
            base = 5*(i-1);
            X_opt(:, base+5) = q(1+lag:T+lag,i);
        end
    
        % Run Policy Test Sims 
        [V_sim, u_sim, p_sim, IVI] = runPolicyTestSims(sysparams, bounds, X_opt, "");
    
        % Total power generated
        p_tot = sum(sum(p_sim)); 
    
        % Store simulation metrics 
        IVI_MC(k, m)     = sum(IVI);
        viol_MC(k,m)     = sum(IVI) > 0;
        P_MC(k, m)       = p_tot;
    end
end 


%% Step 4: HEATMAPS, BENCHMARKING AND REPORTING


% --- choose bin edges (edit as you like) ---
alphaEdges = linspace(0.02,0.40,16);   % 15 bins
DEdges     = linspace(0.25,2.0,16);     % 15 bins

% alphaEdges = linspace(0.02,0.30,16);   % 15 bins
% DEdges     = linspace(0.25, 1.6, 16);     % 15 bins

% --- assign each MC draw to a bin ---
aBin = discretize(amp_MC, alphaEdges);
dBin = discretize(D_MC,   DEdges);

% keep only draws that fall inside bins
ok = ~isnan(aBin) & ~isnan(dBin);
aBin = aBin(ok);
dBin = dBin(ok);

% bin centers (for axes labels)
aCenters = 0.5*(alphaEdges(1:end-1) + alphaEdges(2:end));
dCenters = 0.5*(DEdges(1:end-1)     + DEdges(2:end));

nA = numel(aCenters);
nD = numel(dCenters);

for m = 1:n_pols

    % pull metrics for policy m and restrict to ok draws
    viol = viol_MC(ok, m);    % 0/1
    IVI  = IVI_MC(ok, m);
    Pow  = P_MC(ok, m);

    % count per bin
    Nbin = accumarray([dBin, aBin], 1, [nD, nA], @sum, 0);

    % P(violation) = mean of viol in bin
    Pviol = accumarray([dBin, aBin], viol, [nD, nA], @mean, NaN);

    % E[IVI] = mean of IVI in bin
    EIVI  = accumarray([dBin, aBin], IVI,  [nD, nA], @mean, NaN);

    % E[Power] = mean of power in bin
    EPow  = accumarray([dBin, aBin], Pow,  [nD, nA], @mean, NaN);


    % --- Plot 3 heatmaps for this policy ---
    figure('Name', sprintf('Heatmaps — %s (%s)', policies(m).name, upper(policies(m).code)), ...
           'NumberTitle','off', 'Units','normalized','Position',[0.15 0.15 0.7 0.65]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

    nexttile
    imagesc(aCenters, dCenters, Pviol);
    axis xy; colorbar; grid on;
    xlabel('\alpha (amp1)'); ylabel('D (daysPerEvent)');
    title('P(violation)');

    nexttile
    imagesc(aCenters, dCenters, EIVI);
    axis xy; colorbar; grid on;
    xlabel('\alpha (amp1)'); ylabel('D (daysPerEvent)');
    title('E[IVI]');

    nexttile
    imagesc(aCenters, dCenters, EPow);
    axis xy; colorbar; grid on;
    xlabel('\alpha (amp1)'); ylabel('D (daysPerEvent)');
    title('E[Total Power]');
end

%% ====== FLOW-REGIME SLICES (q0) + HEATMAPS FOR DDU ======

% --- pick the policy index for DDU ---
mDDU = find(strcmp({policies.code}, 'ddu'));
if isempty(mDDU), error('Could not find DDU in policies.code'); end

% --- define LOW / MID / HIGH q0 by quantiles ---
qLo = quantile(q0_MC, 0.25);
qHi = quantile(q0_MC, 0.75);

idxLow  = (q0_MC <= qLo);
idxMid  = (q0_MC >  qLo) & (q0_MC < qHi);
idxHigh = (q0_MC >= qHi);

% helper to compute binned means
binHeat = @(idx, y) localBinnedMean(amp_MC(idx), D_MC(idx), y(idx), alphaEdges, DEdges);

% --- compute heatmaps for each q0 regime (DDU only) ---
regNames = {'Low q0 (<=Q25)','Mid q0 (Q25–Q75)','High q0 (>=Q75)'};
regIdx   = {idxLow, idxMid, idxHigh};

minCount = 10;  % with K=2000 you can afford this; reduces noisy bins

for r = 1:3
    idx = regIdx{r};

    % raw scenario metrics for DDU in this flow regime
    viol = viol_MC(:,mDDU);
    IVI  = IVI_MC(:,mDDU);
    Pow  = P_MC(:,mDDU);

    % binned aggregates over (alpha, D)
    [Pviol, Nbin] = binHeat(idx, viol);  Pviol = Pviol ./ max(Nbin,1);   % mean of 0/1
    [EIVI,  Nbin] = binHeat(idx, IVI);   EIVI  = EIVI  ./ max(Nbin,1);
    [EPow,  Nbin] = binHeat(idx, Pow);   EPow  = EPow  ./ max(Nbin,1);

    % mask sparse bins
    Pviol(Nbin < minCount) = NaN;
    EIVI(Nbin  < minCount) = NaN;
    EPow(Nbin  < minCount) = NaN;

    figure('Name', sprintf('DDU Heatmaps — %s', regNames{r}), ...
        'NumberTitle','off', 'Units','normalized','Position',[0.12 0.18 0.78 0.60]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile
    imagesc(aCenters, dCenters, Pviol); axis xy; colorbar;
    xlabel('\alpha (amp1)'); ylabel('D (days)');
    title(sprintf('P(violation) — %s', regNames{r}));
    grid on;

    nexttile
    imagesc(aCenters, dCenters, EIVI); axis xy; colorbar;
    xlabel('\alpha (amp1)'); ylabel('D (days)');
    title('E[IVI]');
    grid on;

    nexttile
    imagesc(aCenters, dCenters, EPow); axis xy; colorbar;
    xlabel('\alpha (amp1)'); ylabel('D (days)');
    title('E[Total Power]');
    grid on;
end

%% ===== local function (place at bottom of script or in separate file) =====
function [S, N] = localBinnedMean(alpha, D, y, alphaEdges, DEdges)
    aBin = discretize(alpha, alphaEdges);
    dBin = discretize(D,     DEdges);
    ok = ~isnan(aBin) & ~isnan(dBin) & ~isnan(y);

    aBin = aBin(ok); dBin = dBin(ok); y = y(ok);

    nA = numel(alphaEdges)-1;
    nD = numel(DEdges)-1;

    N = accumarray([dBin, aBin], 1, [nD, nA], @sum, 0);
    S = accumarray([dBin, aBin], y, [nD, nA], @sum, 0);
end


%% Difference heatmaps integrated over q0 

% Policy indices
mDIU = find(strcmp({policies.code}, 'diu'));
mDDU = find(strcmp({policies.code}, 'ddu'));
if isempty(mDIU) || isempty(mDDU), error('Could not find DIU/DDU in policies.code'); end

% Assign each draw to a bin
aBin = discretize(amp_MC, alphaEdges);
dBin = discretize(D_MC,   DEdges);
ok   = ~isnan(aBin) & ~isnan(dBin);

aBin = aBin(ok);
dBin = dBin(ok);

% counts per bin
Nbin = accumarray([dBin, aBin], 1, [nD, nA], @sum, 0);

% helper to compute binned means for a metric vector y (length K)
binMean = @(y) accumarray([dBin,aBin], y(ok), [nD,nA], @mean, NaN);

% Binned means (integrated over q0)
Pviol_DIU = binMean(viol_MC(:,mDIU));
Pviol_DDU = binMean(viol_MC(:,mDDU));

EIVI_DIU  = binMean(IVI_MC(:,mDIU));
EIVI_DDU  = binMean(IVI_MC(:,mDDU));

EPow_DIU  = binMean(P_MC(:,mDIU));
EPow_DDU  = binMean(P_MC(:,mDDU));

% Differences (your requested directions)
dPviol = Pviol_DDU - Pviol_DIU;   
dEIVI  = EIVI_DDU - EIVI_DIU;    
dPow   = EPow_DDU  - EPow_DIU;    


% Plot
figure('Name','Difference Heatmaps (q0 integrated)','NumberTitle','off', ...
       'Units','normalized','Position',[0.12 0.18 0.78 0.60]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

nexttile
imagesc(aCenters, dCenters, dPviol); axis xy; colorbar; grid on;
xlabel('\alpha (amp1)'); ylabel('D (days)');
title('DDU - DIU : P(violation)');

nexttile
imagesc(aCenters, dCenters, dEIVI); axis xy; colorbar; grid on;
xlabel('\alpha (amp1)'); ylabel('D (days)');
title('DDU - DIU : E[IVI]');

nexttile
imagesc(aCenters, dCenters, dPow); axis xy; colorbar; grid on;
xlabel('\alpha (amp1)'); ylabel('D (days)');
title('DDU - DIU : E[Total Power]');

fprintf('Simulation complete.\n');
fprintf('Total runtime: %.2f seconds.\n', toc);

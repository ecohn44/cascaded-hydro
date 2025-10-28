%% Author: Eliza Cohn
% Date: October 2025
% Description: Main driver for cascaded hydropower simulations 
% Paper: Optimization of Cascaded Hydroelectric Systems under DDU

clear; clc; close all;
addpath('/Library/gurobi1202/macos_universal2/matlab');

%% ========================================================================
% SECTION 1: DATA LOADING AND PARAMETER DEFINITION
% ========================================================================

% Static Parameters (PowerProd)
eta = .9;          % efficiency of release-energy conversion
rho_w = 1000;      % density of water [kg/m^3]
g = 9.8;           % acceleration due to gravity [m/s^2]
c = eta*rho_w*g/3.6e9; % power prod coefficient
N = 20;             % number of sub-intervals for piecewise linear approx

% Load inflow data
[inflow, modelparams, sysparams] = dataload(N);

%% ========================================================================
% SECTION 2: SIMULATION SETTINGS
% ========================================================================

% Initialize settings (season, linear approximation, uncertainty, bounds)
simSettings = initSimSettings("dry", "pwl", "ddu", "jcc-bon");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));
modelparams.rho = 0.1135; % Calculated offline between (q1_hist, s1 + u1) 

% Date range settings 
D = 2;                        % Simulation duration in days
T = 24*D;                     % Number of simulation hours
lag = 1;                      % Number of lag terms in OLS model
year = 2022;                  % Simulation year

% Compute simulation daterange and inflow series
sim_center_date = datetime(year, 1, 1) + days(modelparams.center_day - 1);
start_date = sim_center_date - hours(T/2) - hours(lag);
end_date   = sim_center_date + hours(T/2 - 1);
inflow_s = inflow(inflow.datetime >= start_date & inflow.datetime <= end_date, :);

% Extract historic inflow timeseries [m3/hr]
% q = [inflow_s.bon_inflow_m3hr];
% q = [parabola_decay(inflow_s.bon_inflow_m3hr(1), T)]';

q0 = inflow_s.bon_inflow_m3hr(1);      % your baseline level
% q = makeInflowPulse(q0, T, lag, [0.2 0.6], 0.20, 0.80, 2, 2);
q = q0*ones(T+1,1);

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% ========================================================================

[model, obj, X, std_hat, phi_vals, alpha_vals, U_eff] = optimization(T, N, c, q, lag, ...
    simSettings.framework, simSettings.bounds, modelparams, sysparams);

% Extract q2 reference inflow
q(:,2) = [0; X(:,3) + X(:,4)];

% Calculate percent difference (step change) for flow bounds (flow_max 02)
pct_diff = 100 * diff(U_eff(:,3)) ./ U_eff(1:end-1,3);

% Calculate percent difference (from max) for flow bounds (flow_max 02)
max_V2 = sysparams(2).max_V;
pct_diff_from_max = 100 * (max_V2 - U_eff(:,3)) ./ max_V2;

%% ========================================================================
% SECTION 4: PLOTTING
% ========================================================================

% Toggle for creating folder and plotting
make_dir = false; % Set to true to enable directory creation and plotting
printplot = false; 

% Make plot directory for current simulation run 
if make_dir
    dir_path = "./plots/";
    stamp = datestr(now,'mm-dd-yyyy HH.MM.SS');
    path = fullfile(dir_path, stamp + " " + simSettings.season + " " ...
        + simSettings.framework + ...
        " T=" + string(T));
    mkdir(path)
end

% Plot simulation behavior for all units
simPlots(path, X, U_eff, q, sysparams, T, c, lag, printplot);

% Make simulation result storage folder
results_dir = "./results/";
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Store simulation results 
for i = 1:numel(sysparams)
    sp = sysparams(i);
    fname = sprintf('results_unit%d_%s.mat', sp.unit, lower(simSettings.framework));
    save(fullfile(results_dir, fname), ...
        'X', 'U_eff', 'std_hat', 'q', 'sysparams', 'T', 'c', 'lag', '-v7');
end

%% ========================================================================
% SECTION 5: MONTE CARLO SIMS
% ========================================================================

monte_carlo = false;

if monte_carlo
    fprintf('Running Monte Carlo Sims.\n');
    [V1, V2] = runMonteCarloSims(sysparams, simSettings.bounds, std_hat, X, path, printplot);
end

fprintf('Simulation complete.\n');

% Show dynamic slack allocation in chance constraints 
if simSettings.bounds == "jcc-ssh"
    figure;
    plot(alpha_vals(:,1)*100, 'b-', 'LineWidth',1.5); hold on;
    plot(alpha_vals(:,2)*100, 'r--', 'LineWidth',1.5);
    yline(100*(1/2), 'k:', 'Bonferroni 50/50');
    xlabel('Time step');
    ylabel('Slack allocation (%)');
    legend('Reservoir 1','Reservoir 2','Location','best');
    title('Adaptive Risk Allocation (SSH vs. Bonferroni)');
    grid on;
end 


function q = makeInflowPulse(q0, T, lag, t0, amp1, amp2, w1, w2)
% makeInflowPulse  Two positive pulses in q1 at arbitrary positions.
% q0   : scalar baseline (m^3/hr) OR vector of length T+lag
% T,lag: horizon & lag (output length = T+lag)
% t0   : pulse positions:
%        - if scalar <1: one pulse at fraction t0 of horizon
%        - if scalar >=1: one pulse at absolute index
%        - if 2-vector: two pulses; each entry may be fraction (<1) or index
% amp1 : fractional bump for pulse #1 (e.g., 0.2 => +20%)
% amp2 : fractional bump for pulse #2 (e.g., 0.8 => +80%)
% w1,w2: widths (hours) for pulses #1 and #2
%
% Returns: q (T+lag x 1) nonnegative

    % --- Baseline ---
    if isscalar(q0)
        q = q0 * ones(T+lag,1);
    else
        assert(numel(q0) == T+lag, 'q0 must be scalar or length T+lag');
        q = q0(:);
    end

    % --- Normalize positions to absolute indices in 1..T+lag ---
    if nargin < 4 || isempty(t0), t0 = [0.2 0.8]; end     % default 20% and 80%
    t0 = t0(:)';                                           % row
    if numel(t0) == 1, t0 = [t0]; end

    toIndex = @(x) ( x < 1 ) .* (lag + round(x*T)) + ( x >= 1 ) .* round(x);
    t_idx = arrayfun(@(x) max(1, min(T+lag, toIndex(x))), t0);

    % --- Widths (allow 0 => no pulse) ---
    w1 = max(0, round(w1));  w2 = max(0, round(w2));
    amp1 = max(0, amp1);     amp2 = max(0, amp2);

    % --- Helper: apply one centered box pulse (multiplicative) ---
    function applyPulse(center_idx, width, amp)
        if width <= 0 || amp <= 0, return; end
        half = floor((width-1)/2);
        i_start = max(1, center_idx - half);
        i_end   = min(T+lag, i_start + width - 1);
        q(i_start:i_end) = q(i_start:i_end) .* (1 + amp);
    end

    % Pulse #1 at t_idx(1), Pulse #2 at t_idx(2) if present
    applyPulse(t_idx(1), w1, amp1);
    if numel(t_idx) >= 2
        applyPulse(t_idx(2), w2, amp2);
    end

    % --- Nonnegative safeguard ---
    q = max(q, 0);
end

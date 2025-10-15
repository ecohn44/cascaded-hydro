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
simSettings = initSimSettings("dry", "pwl", "diu", "icc");

% Extract forecasting coefficients 
modelparams = modelparams(strcmp({modelparams.season}, simSettings.season));
modelparams.rho = 0.1135; % Calculated offline between (q1_hist, s1 + u1) 

% Date range settings 
D = 5;                        % Simulation duration in days
% T = 24*D;                   % Number of simulation hours
T = 48;
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
t0   = round(0.5*T) + lag;               % center pulse mid-horizon
q = makeInflowPulse(q0, T, lag, 0.35, 0.6, 0.4, 2, 2);

fprintf('Data loading complete.\n');

%% ========================================================================
% SECTION 3: OPTIMIZATION FRAMEWORK
% ========================================================================

[model, obj, X, std_hat, phi_vals, alpha_vals, U_eff] = optimization(T, N, c, q, lag, ...
    simSettings.framework, simSettings.bounds, modelparams, sysparams);

% Extract q2 reference inflow
q(:,2) = [0; X(:,3) + X(:,4)];


%% ========================================================================
% SECTION 4: PLOTTING
% ========================================================================

% Toggle for creating folder and plotting
make_dir = true; % Set to true to enable directory creation and plotting
printplot = true; 

if make_dir
    dir_path = "./plots/";
    stamp = datestr(now,'mm-dd-yyyy HH.MM.SS');
    path = fullfile(dir_path, stamp + " " + simSettings.season + " " ...
        + simSettings.framework + ...
        " T=" + string(T));
    mkdir(path)
end

% Plot simulation behavior for all units
simPlots(path, X, q, sysparams, T, c, lag, printplot);

% Save results
save(sprintf('unit2_%s.mat', lower(simSettings.framework)), 'X','U_eff','std_hat','sysparams','-v7');

%% ========================================================================
% SECTION 5: MONTE CARLO SIMS
% ========================================================================

fprintf('Running Monte Carlo Sims.\n');

[V1, V2] = runMonteCarloSims(sysparams, simSettings.bounds, std_hat, X, path, printplot);

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


function q = makeInflowPulse(q0, T, lag, t0, amp_up, amp_dn, w_up, w_dn)
% makeInflowPulse  Create an upstream inflow pulse (q1) mid-horizon.
% q0      : scalar baseline (m^3/hr) OR vector of length T+lag
% T, lag  : horizon length and model lag (q must be length T+lag)
% t0      : pulse location. If t0<1, interpreted as fraction of horizon (0..1);
%           otherwise treated as absolute index in 1..T+lag.
% amp_up  : fractional bump during expansion (e.g., 0.6 => +60%)
% amp_dn  : fractional dip during contraction (e.g., 0.4 => −40%)
% w_up    : width (hours) of expansion segment
% w_dn    : width (hours) of contraction segment (immediately after)
%
% Returns:
%   q     : vector length T+lag with the pulse applied (nonnegative)

    % Baseline
    if isscalar(q0)
        q = q0 * ones(T+lag,1);
    else
        assert(numel(q0) == T+lag, 'q0 must be scalar or length T+lag');
        q = q0(:);
    end

    % --- Interpret t0 ---
    if t0 < 1                   % allow fractional placement in [0,1]
        t0 = lag + round(t0 * T);
    end
    t0   = max(1, min(T+lag, round(t0)));
    w_up = max(0, round(w_up));
    w_dn = max(0, round(w_dn));

    % Expansion (centered at t0)
    if w_up > 0
        half = floor((w_up-1)/2);
        i_up_start = max(1, t0 - half);
        i_up_end   = min(T+lag, i_up_start + w_up - 1);
        q(i_up_start:i_up_end) = q(i_up_start:i_up_end) .* (1 + amp_up);
        last_up = i_up_end;
    else
        last_up = t0;
    end

    % Contraction (right after expansion)
    if w_dn > 0
        i_dn_start = min(T+lag, last_up + 1);
        i_dn_end   = min(T+lag, i_dn_start + w_dn - 1);
        if i_dn_end >= i_dn_start
            q(i_dn_start:i_dn_end) = q(i_dn_start:i_dn_end) .* (1 - amp_dn);
        end
    end

    % Nonnegative safeguard
    q = max(q, 0);
end

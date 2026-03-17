%% ========================================================================
% Author: Eliza Cohn
% Date: March 2026
% Description: Monte Carlo violation test (aligned with phi_sweep)
% =========================================================================

clc; clear; close all;

K = 10000;         % number of MC simulations
std = 0.015;      % noise std
rng(0)

% Load precomputed results
load('./mats/drought_sweep_jcc-bon.mat');
phi_bon = phi_sweep; X_bon = X_sweep; q = q_sweep;

load('./mats/drought_sweep_jcc-ssh.mat');
phi_ssh = phi_sweep; X_ssh = X_sweep;

n   = 3;
N   = 40;        % piecewise intervals
[modelparams, sysparams, seasonparams] = dataload(n, N);

% Simulation horizon
D = 3.5;         
T = 24*D;        
lag = 1;

Vmin = sysparams(1).min_V;   % minimum storage
Vmax = sysparams(1).max_V;   % maximum storage

alpha_set = 0.1:0.05:0.4;
baseline_set = 0.02:0.0025:0.045;

empirical_phi_ssh = zeros(length(alpha_set), length(baseline_set));

for a = 1:length(alpha_set)
    for d = 1:length(baseline_set)
        
        X_opt = X_ssh{a,d};
        q_cell = q{a,d};

        min_success_per_run = zeros(K,1);

        for k = 1:K
            % Sample the full trajectory (T x n) for this run
            noise = std * randn(size(q_cell));   % T x n
            q_sample = q_cell + noise;
            q_sample = max(q_sample, 0);         % avoid negative inflows

            % Overwrite inflows in X_opt for the full trajectory
            X_run = X_opt;
            for i = 1:n
                base = 5*(i-1);
                X_run(:, base+5) = q_sample(:, i);
            end

            % Run full trajectory simulation
            [V_sim, ~, ~, ~] = runPolicyTestSims(sysparams, 'jcc-ssh', X_run, "");

            % Check violations at each timestep
            success_t = all(V_sim >= Vmin & V_sim <= Vmax, 2);  % T x 1

            % Minimum over timesteps for this run
            min_success_per_run(k) = min(success_t);
        end

        % Empirical phi: fraction of runs that satisfy all timesteps
        empirical_phi_ssh(a,d) = mean(min_success_per_run);
    end
end

% Compare trends
disp(phi_ssh)
disp(empirical_phi_ssh)

phi_diff = phi_ssh - empirical_phi_ssh;
phi_diff_pct = 100 * (phi_ssh - empirical_phi_ssh) ./ phi_ssh;

figure;
imagesc(baseline_set, alpha_set, phi_diff_pct);
xlabel('Baseline Flow (q_0)');
ylabel('Drought Intensity (\alpha)');
title('Predicted SSH \phi - Empirical \phi');

% Custom green-to-red colormap
axis xy; grid on;
cb = colorbar;
cb.Label.String = 'Percent Difference (%)';
cb.Label.FontWeight = 'bold';
colormap(parula);       


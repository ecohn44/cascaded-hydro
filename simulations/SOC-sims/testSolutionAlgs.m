%% ========================================================================
% Author: Eliza Cohn
% Date: March 2026
% Description: Comparing the performance of BON to SSH on drought scenarios
% =========================================================================
clc; clear; close all; 

load('./mats/drought_sweep_jcc-bon.mat');
phi_bon = phi_sweep; X_bon = X_sweep; q = q_sweep;

load('./mats/drought_sweep_jcc-ssh.mat');
phi_ssh = phi_sweep; X_ssh = X_sweep;

n   = 3;
eps = 0.05;
[~, sysparams, ~] = dataload(n, 40);

% Storage
IVI_bon    = nan(length(alpha_set), length(baseline_set));
IVI_ssh    = nan(length(alpha_set), length(baseline_set));
energy_bon = nan(length(alpha_set), length(baseline_set));
energy_ssh = nan(length(alpha_set), length(baseline_set));


for a = 1:length(alpha_set)
    for d = 1:length(baseline_set)

        q_cell = q{a, d};

        % BON: use own policy if converged, else nearest neighbor
        if ~isnan(phi_bon(a,d))
            X_bon_cell = X_bon{a, d};
        else
            [a_nb, d_nb] = nearestNeighbor(phi_bon, a, d);
            X_bon_cell   = X_bon{a_nb, d_nb};
        end

        for i = 1:n
            base = 5*(i-1);
            X_bon_cell(:, base+5) = q_cell(:, i);
        end

        [V_bon, u_bon, p_bon, IVI_bon_val] = runPolicyTestSims(sysparams, 'jcc-bon', X_bon_cell, "");
        IVI_bon(a,d)    = sum(rescale(IVI_bon_val, "volume"));
        h_phys_bon = rescale(V_bon, "head");
        u_phys_bon = rescale(u_bon, "release");
        energy_bon(a,d) = sum(rescale(h_phys_bon.*u_phys_bon, "power"), 'all');

        % SSH: always own policy
        X_ssh_cell = X_ssh{a, d};
        for i = 1:n
            base = 5*(i-1);
            X_ssh_cell(:, base+5) = q_cell(:, i);
        end
        
        [V_ssh, u_ssh, p_ssh, IVI_ssh_val] = runPolicyTestSims(sysparams, 'jcc-ssh', X_ssh_cell, "");
        IVI_ssh(a,d)    = sum(rescale(IVI_ssh_val, "volume"));
        h_phys_ssh = rescale(V_ssh, "head");
        u_phys_ssh = rescale(u_ssh, "release");
        energy_ssh(a,d) = sum(rescale(h_phys_ssh.*u_phys_ssh, "power"), 'all');

    end
end

% Plot
plotPolicyComparison(IVI_bon, IVI_ssh, energy_bon, energy_ssh, alpha_set, baseline_set);

IVI_diff    = IVI_bon    - IVI_ssh;
energy_diff = energy_ssh - energy_bon;

plotHeat(alpha_set, baseline_set, IVI_diff', 0, 0, "IVI_diff", ...
    'Drought Intensity (\alpha)', 'Baseline Flow (q_0)');

plotHeat(alpha_set, baseline_set, energy_diff', 0, 0, "energy", ...
    'Drought Intensity (\alpha)', 'Baseline Flow (q_0)');



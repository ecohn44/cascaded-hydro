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


%% 2D Line Plot
figure; hold on; grid on;

energy_ssh = energy_ssh/1e3;
energy_bon = energy_bon/1e3;

% Partition low and high flows 
nq = size(energy_ssh,2);
low_idx  = 1:floor(nq/2);
high_idx = floor(nq/2)+1:nq;

c_bon = [0 0.4470 0.7410];      % benchmark (blue)
c_bon_dark = [0 0.2 0.5];
c_ssh = [0.2 0.8 0.2];    % light green shading
c_ssh_dark  = [0 0.5 0];        % dark green mean line

mean_ssh_low  = mean(energy_ssh(:,low_idx), 2, 'omitnan');
mean_ssh_high = mean(energy_ssh(:,high_idx),2, 'omitnan');
mean_bon_low  = mean(energy_bon(:,low_idx), 2, 'omitnan');
mean_bon_high = mean(energy_bon(:,high_idx),2, 'omitnan');


std_ssh_low  = std(energy_ssh(:,low_idx), 0, 2);
std_ssh_high = std(energy_ssh(:,high_idx),0, 2);
std_bon_low  = std(energy_bon(:,low_idx), 0, 2);
std_bon_high = std(energy_bon(:,high_idx),0, 2);


% BON Low Flow
fill([alpha_set, fliplr(alpha_set)], ...
     [mean_bon_low'+std_bon_low', fliplr(mean_bon_low'-std_bon_low')], ...
     c_bon, 'FaceAlpha',0.2, 'EdgeColor','none');

% BON High Flow
fill([alpha_set, fliplr(alpha_set)], ...
     [mean_bon_high'+std_bon_high', fliplr(mean_bon_high'-std_bon_high')], ...
     c_bon, 'FaceAlpha',0.1, 'EdgeColor','none');  % lighter for distinction

% SSH Low Flow
fill([alpha_set, fliplr(alpha_set)], ...
     [mean_ssh_low'+std_ssh_low', fliplr(mean_ssh_low'-std_ssh_low')], ...
     c_ssh, 'FaceAlpha',0.2, 'EdgeColor','none');

% SSH High Flow
fill([alpha_set, fliplr(alpha_set)], ...
     [mean_ssh_high'+std_ssh_high', fliplr(mean_ssh_high'-std_ssh_high')], ...
     c_ssh, 'FaceAlpha',0.1, 'EdgeColor','none');

% Mean lines
plot(alpha_set, mean_bon_low,  '-o','LineWidth',2.5,'Color',c_bon_dark, 'MarkerFaceColor',c_bon_dark,'MarkerSize',10)
plot(alpha_set, mean_bon_high, '-s','LineWidth',2.5,'Color',c_bon_dark, 'MarkerSize',10)
plot(alpha_set, mean_ssh_low,  '-o','LineWidth',2.5,'Color',c_ssh_dark, 'MarkerFaceColor',c_ssh_dark,'MarkerSize',10)
plot(alpha_set, mean_ssh_high, '-s','LineWidth',2.5,'Color',c_ssh_dark, 'MarkerSize',10)

ax = gca;                  % get current axes
ax.FontSize = 14;           % increase tick labels (x and y)
ax.LineWidth = 1.5;         % optional: thicker axes lines for clarity

xlabel('Drought Intensity \alpha', FontSize=16)
ylabel('Total Energy (GWh)', FontSize=16)
xlim([min(alpha_set) max(alpha_set)])

% Legend
hMean_ssh = [plot(NaN,NaN,'-o','Color',c_ssh_dark,'MarkerFaceColor',c_ssh_dark,'LineWidth',2), ...
         plot(NaN,NaN,'-s','Color',c_ssh_dark,'LineWidth',2)];
hBand_ssh = [fill(NaN,NaN,c_ssh,'FaceAlpha',0.2,'EdgeColor','none')];

hMean_bon = [plot(NaN,NaN,'-o','Color',c_bon_dark,'MarkerFaceColor',c_bon_dark,'LineWidth',2), ...
         plot(NaN,NaN,'-s','Color',c_bon_dark,'LineWidth',2)];

hBand_bon = [fill(NaN,NaN,c_bon,'FaceAlpha',0.2,'EdgeColor','none')];

hLegend = legend([hMean_ssh hBand_ssh hMean_bon hBand_bon], {'SSH Low Flow','SSH High Flow', 'SSH ±1σ'...
                       'BON Low Flow','BON High Flow','BON ±1σ'});

% Make it larger
hLegend.FontSize = 16;          % increase font size
hLegend.Box = 'off';            % optional, remove box

% Place below figure
hLegend.Orientation = 'horizontal';
hLegend.Location = 'southoutside';  % below the axes
hLegend.NumColumns = 3;             % 2 columns: SSH vs BON

grid off

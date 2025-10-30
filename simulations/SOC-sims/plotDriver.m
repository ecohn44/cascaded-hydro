% clear; clc; close all;

% --- Load results from the universal results folder ---
S_diu = load('./results/results_unit2_diu.mat');
S_ddu = load('./results/results_unit2_ddu.mat');

% --- Plot options ---
opts = struct('showBands', true, 'showDelta', false, 'titleTag', 'dry');

% ------------------------------------------------------------
% DIU vs DDU overlay for Unit 2 (SoC + sigma) using V_eff
% assumes you already ran:
%   S_diu = load('./results/results_unit2_diu.mat');
%   S_ddu = load('./results/results_unit2_ddu.mat');
% and you have: S_*.X, S_*.std_hat, S_*.V_eff, S_ddu.sysparams
% ------------------------------------------------------------

X_diu    = S_diu.X;
X_ddu    = S_ddu.X;
std_diu  = S_diu.std_hat;
std_ddu  = S_ddu.std_hat;
Veff_diu = S_diu.V_eff;
Veff_ddu = S_ddu.V_eff;
s        = S_ddu.sysparams;     % sysparams array
savePath = './plots';
opts     = struct('titleTag','dry');

if ~exist(savePath,'dir')
    mkdir(savePath);
end

T  = size(X_ddu,1);
tt = (1:T)';

% ------------------- unit 2 series -------------------
V2_diu = X_diu(:,6);
V2_ddu = X_ddu(:,6);
q2_diu = X_diu(:,10); %#ok<NASGU>
q2_ddu = X_ddu(:,10); %#ok<NASGU>

Vmin2  = s(2).min_V;
Vmax2  = s(2).max_V;
span2  = Vmax2 - Vmin2;

soc_diu = (V2_diu - Vmin2)/span2;
soc_ddu = (V2_ddu - Vmin2)/span2;

% event = time of max downstream std in DDU
[~, t_evt] = max(std_ddu(:,2));

% =====================================================
% FIGURE 1: sigma_2(t)
% =====================================================
f1 = figure('Position',[100 100 850 350]); hold on; grid on;
plot(tt, std_diu(:,2), 'LineWidth',1.5, 'DisplayName','DIU \sigma_2');
plot(tt, std_ddu(:,2), 'LineWidth',1.5, 'DisplayName','DDU \sigma_2');
xline(t_evt, 'k--', 'spike', 'LabelVerticalAlignment','bottom');
xlabel('time step');
ylabel('\sigma_{q2} (normalized)');
title(['Unit 2 inflow uncertainty (DIU vs DDU) ' opts.titleTag]);
legend('Location','best');

saveas(f1, fullfile(savePath, 'u2_sigma_diu_vs_ddu.png'));

% =====================================================
% FIGURE 2: SoC(t) + corridors from V_eff
% V_eff format is: [V1_max V1_min V2_max V2_min]
% =====================================================
V2max_diu = Veff_diu(:,3);
V2min_diu = Veff_diu(:,4);
V2max_ddu = Veff_ddu(:,3);
V2min_ddu = Veff_ddu(:,4);

soc2max_diu = (V2max_diu - Vmin2)/span2;
soc2min_diu = (V2min_diu - Vmin2)/span2;
soc2max_ddu = (V2max_ddu - Vmin2)/span2;
soc2min_ddu = (V2min_ddu - Vmin2)/span2;

f2 = figure('Position',[100 480 1100 480]); hold on; grid on;
yline(0,'k:');
yline(1,'k:');

% DIU corridor (green-ish)
patch([tt; flipud(tt)], [soc2min_diu; flipud(soc2max_diu)], [0.75 1.0 0.75], ...
      'EdgeColor','none', 'FaceAlpha',0.35, 'DisplayName','DIU corridor');

% DDU corridor (red-ish)
patch([tt; flipud(tt)], [soc2min_ddu; flipud(soc2max_ddu)], [1.0 0.7 0.6], ...
      'EdgeColor','none', 'FaceAlpha',0.35, 'DisplayName','DDU corridor');

% actual SoC lines
plot(tt, soc_diu, 'Color',[0.1 0.5 0.1], 'LineWidth',2, 'DisplayName','DIU SoC');
plot(tt, soc_ddu, 'Color',[0.7 0.0 0.0], 'LineWidth',2, 'DisplayName','DDU SoC');

xline(t_evt, 'k--', 'risk event', 'LabelVerticalAlignment','bottom');

xlabel('time step');
ylabel('SoC_2 = (V_2 - V_{min})/(V_{max} - V_{min})');
title(['Unit 2 SoC: DIU vs DDU ' opts.titleTag]);
legend('Location','best');

saveas(f2, fullfile(savePath, 'u2_soc_diu_vs_ddu.png'));

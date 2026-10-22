%MAKE_PRICE_FIGURE  The settlement price signal: BPA PF-24 Load Shaping rate.
%
%   Run it:  >> make_price_figure        (script, not a function)
%
%   Writes TWO figures, sized and lettered for the paper:
%     figures_out/bpa_price_signal.{pdf,png}
%         top    -- the hourly schedule over WY2024
%         bottom -- a week in June (freshet) and a week in August (peak)
%     figures_out/bpa_price_tariff.{pdf,png}
%         the published monthly tariff the schedule is expanded from
%
%   No titles: label the panels with \subfloat in the manuscript so the
%   caption carries the text.
%
%   Source: Bonneville Power Administration, 2024 Power Rate Schedules and
%   General Rate Schedule Provisions (BP-24), Schedule PF-24 Priority Firm
%   Power, sec. 2.1.3.1 "Load Shaping Rate", p. 8. Rates in mills/kWh, which
%   are numerically equal to $/MWh, and are used unmodified.
%
%   Heavy Load Hours: hours ending 0700-2200 Pacific Prevailing Time,
%   Monday-Saturday, excluding NERC holidays. All other hours, all of Sunday
%   and all NERC holidays are Light Load Hours.

clear; clc;

%% ---------------------------------------------------------------- CONFIG
here = fileparts(mfilename('fullpath'));
ROOT = fileparts(here);                        % the repository root
DATA = fullfile(ROOT, 'data');
FIGD = fullfile(ROOT, 'figures_out');
if ~isfolder(FIGD), mkdir(FIGD); end

HLH_C = [  0 114 178]/255;                     % validated CVD-safe pair
LLH_C = [  0 158 115]/255;
FS    = 16;                                    % print size

WK_JUN = datetime(2024,6,3);                   % Monday
WK_AUG = datetime(2024,8,5);                   % Monday

%% ------------------------------------------------------------ LOAD DATA
P = readtable(fullfile(DATA,'bpa_pf24_lmp_wy2024.csv'), 'VariableNamingRule','preserve');
P.dt = datetime(P.datetime, 'InputFormat','yyyy-MM-dd HH:mm');
R = readtable(fullfile(DATA,'bpa_pf24_rates.csv'), 'VariableNamingRule','preserve');

wyOrder = [10 11 12 1 2 3 4 5 6 7 8 9];
lab = {'Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'};
hlh = zeros(12,1); llh = zeros(12,1);
for k = 1:12
    r = R(R.month == wyOrder(k), :);
    hlh(k) = r.hlh_usd_per_mwh;  llh(k) = r.llh_usd_per_mwh;
end
annMean = mean(P.price_usd_per_mwh);

%% ============================ FIGURE 1: schedule + two example weeks ====
f1 = figure('Units','inches','Position',[1 1 7.16 2.6],'Color','w');
ax = axes(f1); hold(ax,'on'); grid(ax,'off');

stairs(ax, P.dt, P.price_usd_per_mwh, 'Color', HLH_C, 'LineWidth', 0.8);
ylim(ax,[0 80]); xlim(ax, [min(P.dt) max(P.dt)]);
xticks(ax, datetime(2023,10,1) + calmonths(0:11));
xtickformat(ax, 'MMM');
ax.XAxis.SecondaryLabel.Visible = 'off';       % drop the "2023-2024" row
ylabel(ax,'$/MWh');
styleAxes(ax, FS);

exportgraphics(f1, fullfile(FIGD,'bpa_price_signal.pdf'), 'ContentType','vector');
exportgraphics(f1, fullfile(FIGD,'bpa_price_signal.png'), 'Resolution', 400);
%% ============================ FIGURE 2: the published tariff ===========
f2 = figure('Units','inches','Position',[1 1 7.16 2.6],'Color','w');
ax = axes(f2); hold(ax,'on'); grid(ax,'off');
b = bar(ax, 1:12, [hlh llh], 0.82, 'EdgeColor','none');
b(1).FaceColor = HLH_C;  b(1).DisplayName = 'HLH';
b(2).FaceColor = LLH_C;  b(2).DisplayName = 'LLH';
yline(ax, annMean, ':', 'Color',[0.45 0.45 0.45], 'LineWidth',1.5, ...
      'DisplayName', 'Annual Mean');
xticks(ax,1:12); xticklabels(ax, lab); xlim(ax,[0.4 12.6]);
ylabel(ax,'$/MWh');
legend(ax, 'Location','northwest', 'Box','off', 'NumColumns',3);
styleAxes(ax, FS);

exportgraphics(f2, fullfile(FIGD,'bpa_price_tariff.pdf'), 'ContentType','vector');
exportgraphics(f2, fullfile(FIGD,'bpa_price_tariff.png'), 'Resolution', 400);

%% ---------------------------------------------------------------- REPORT
fprintf('annual mean $%.2f | min $%.2f | max $%.2f | HLH %.1f%% of hours\n', ...
    annMean, min(P.price_usd_per_mwh), max(P.price_usd_per_mwh), ...
    100*mean(P.is_hlh));
fprintf('seasonal spread %.2fx  (Aug HLH $%.2f / Jun LLH $%.2f)\n', ...
    max(hlh)/min(llh), max(hlh), min(llh));
fprintf('wrote bpa_price_signal.{pdf,png} and bpa_price_tariff.{pdf,png} to %s\n', FIGD);

%% ------------------------------------------------------------- helpers
function styleAxes(ax, FS)
ax.FontSize = FS; ax.FontName = 'Times New Roman';
ax.LineWidth = 0.7; ax.XColor = [0.25 0.25 0.25]; ax.YColor = [0.25 0.25 0.25];
ax.GridColor = [0.84 0.84 0.83]; ax.GridAlpha = 1; ax.Box = 'on';
ax.XLabel.FontSize = FS; ax.YLabel.FontSize = FS;
end

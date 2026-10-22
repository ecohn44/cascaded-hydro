%MAKE_TRACKING_FIGURE  Monthly mean storage against the seasonal reference,
%                      for each tracking weight theta.
%
%   Run it:  >> make_tracking_figure        (script, not a function)
%
%   Replaces the archived Python version
%       archive/2026-09-25_cleanup/paper_figures_before_move/
%           make_traj_tracking_symmetric.py
%   which predates the BPA price signal and the WY2024 grid. Everything is
%   read from the finished BPA_HL2 runs; nothing is refitted or re-solved.
%
%   Envelope and reference come from load_case, not from a saved pickle, so
%   they are guaranteed to be the same ones the controller was given.
%
%   Output: figures_out/traj_tracking_theta.{pdf,png}

clear; clc;

%% ---------------------------------------------------------------- CONFIG
here = fileparts(mfilename('fullpath'));
ROOT = fileparts(here);                        % the repository root
addpath(ROOT); addpath(fullfile(ROOT,'functions'));
FIGD = fullfile(ROOT,'figures_out');
if ~isfolder(FIGD), mkdir(FIGD); end

THETAS = [100 200 400 800];        % the headline ladder; see note below
JOB    = 's01_e010_DDU_k100_BON.mat';   % scenario 1, eps = 0.01, DDU, BON
UNIT   = 'JDA';                    % 'MCN' | 'JDA' | 'TDA' | 'BON'
YEAR   = 2024;  FLOORPCT = 0.15;  HALFWIN = 7;
T0     = datetime(2023,10,1);      % WY2024 starts here

% House conventions, matching window_figure_TDA.m and make_fig_tracking_demo.m:
% the operating floor is a grey fill under a black dashed outline, the
% reference takes the operating-envelope blue, and theta -- an ordered
% magnitude -- takes a single-hue sequential green ramp, light to dark.
FILL_C  = [0.90 0.89 0.87];        % P20 floor fill
FLOOR_C = [0 0 0];                 % P20 floor outline, black dashed
REF_C   = [ 42 120 214]/255;       % envelope blue
RAMP = [162 208 124; 108 174 72; 67 133 44; 15 58 11]/255;
FS   = 11;

%% ------------------------------------------------- envelope and reference
C = load_case(ROOT, YEAR, FLOORPCT, HALFWIN);
iu = find(strcmp(C.units, UNIT));
assert(~isempty(iu), 'UNIT must be one of %s', strjoin(C.units,', '));

%% ---------------------------------------------------------- trajectories
V = cell(numel(THETAS),1);
for k = 1:numel(THETAS)
    f = fullfile(ROOT,'results','runs', ...
                 sprintf('BPA_HL2_TH%07d',THETAS(k)),'jobs',JOB);
    assert(isfile(f), 'missing run: %s', f);
    S = load(f,'small');
    V{k} = S.small.V;
end
T = size(V{1},1);

%% --------------------------------------------------- monthly means, Jan-Dec
% The record is one water year (Oct 2023 - Sep 2024), so Oct-Dec precede
% Jan-Sep chronologically; the axis is ordered by calendar month.
MCM = 1.23348;                     % 10^6 m^3 per kaf (1 af = 1233.48 m^3)
dt  = T0 + hours(0:T-1);
mo  = month(dt(:));
cal = 1:12;
MON = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
mm  = @(v) MCM * arrayfun(@(m) mean(v(mo == m), 'omitnan'), cal(:));

lo  = mm(C.qLow(1:T,iu));
rf  = mm(C.ref(1:T,iu));
mV  = cellfun(@(v) mm(v(:,iu)), V, 'UniformOutput', false);

%% ---------------------------------------------------------------- FIGURE
f = figure('Units','inches','Position',[1 1 6.2 3.3],'Color','w');
ax = axes(f); hold(ax,'on'); grid(ax,'off');
x = 1:12;

% shade the region below the floor
yl = [floor(min(lo)/5)*5 - 5, ceil(max(cellfun(@max,mV))/5)*5 + 5];
yl(1) = max(yl(1), 0);
patch(ax, [x fliplr(x)], [lo(:)' repmat(yl(1),1,12)], FILL_C, ...
      'EdgeColor','none', 'HandleVisibility','off');

plot(ax, x, lo, '--', 'Color',FLOOR_C, 'LineWidth',1.1, 'DisplayName','P10 floor');
plot(ax, x, rf, '-', 'Color',REF_C, 'LineWidth',2.0, 'DisplayName','Reference');
for k = 1:numel(THETAS)
    plot(ax, x, mV{k}, '-o', 'Color',RAMP(k,:), 'LineWidth',1.5, ...
         'MarkerSize',5, 'MarkerFaceColor',RAMP(k,:), 'MarkerEdgeColor',RAMP(k,:), ...
         'DisplayName', sprintf('\\theta = %d', THETAS(k)));
end

xticks(ax,x); xticklabels(ax,MON); xtickangle(ax,-30); xlim(ax,[1, 12]); ylim(ax,yl);
ylabel(ax,'Monthly Mean Storage (10^6 m^3)');
legend(ax,'Location','southoutside','Box','off','NumColumns',3);
ax.FontSize = FS; ax.FontName = 'Times New Roman';
ax.LineWidth = 0.7; ax.XColor = [0.25 0.25 0.25]; ax.YColor = [0.25 0.25 0.25];
ax.GridColor = [0.86 0.86 0.85]; ax.GridAlpha = 1; ax.Box = 'on';

exportgraphics(f, fullfile(FIGD,'traj_tracking_theta.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(FIGD,'traj_tracking_theta.png'), 'Resolution', 400);

%% ---------------------------------------------------------------- REPORT
fprintf('%s, %s   (storage in 10^6 m^3)\n', UNIT, strrep(JOB,'_','\_'));
fprintf('  %-16s mean %6.2f  min %6.2f  max %6.2f\n', 'P20 floor', mean(lo), min(lo), max(lo));
fprintf('  %-16s mean %6.2f  min %6.2f  max %6.2f\n', 'reference', mean(rf), min(rf), max(rf));
for k = 1:numel(THETAS)
    v = mV{k};
    fprintf('  theta %-10d mean %6.2f  min %6.2f  max %6.2f   |dev| %5.2f\n', ...
            THETAS(k), mean(v), min(v), max(v), mean(abs(v - rf)));
end
fprintf('wrote %s\n', fullfile(FIGD,'traj_tracking_theta.pdf'));

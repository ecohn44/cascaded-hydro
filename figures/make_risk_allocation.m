%MAKE_RISK_ALLOCATION  How SSH spends the chance-constraint budget, DDU only.
%
%   Run it:  >> make_risk_allocation        (script, not a function)
%
%   One panel. The horizontal axis is how much room the tightest unit has
%   above its P20 floor; the vertical axis is that unit's share of the joint
%   risk budget. Two operating regimes fall out of the data rather than being
%   imposed:
%
%     even       nothing is near a floor, every share sits on eps/n, and the
%                allocation carries no information
%     dynamic    a floor is in reach, the shares fan out, and in practice one
%                project carries almost the whole budget
%
%   The even branch is explicit in realtimeSimulation.m: when the joint
%   reliability has no measurable sensitivity to storage,
%       sum(abs(gV1)) <= 1e-14  ->  alpha = ones(1,n)/n
%   so eps/n is a principled fallback taken exactly when there is no risk to
%   allocate, NOT a limit the gradient ratio converges to. Say it that way in
%   the caption; a reviewer who opens the code will check.
%
%   Bonferroni is not simulated here. Its allocation IS eps/n by construction,
%   which is the dashed line, so the figure compares the two mechanisms
%   directly without needing the BON runs.
%
%   Source: the SSH trajectories saved at SaveTrajEpsilon = 0.05. Those are
%   the only epsilon with per-hour risk_allocation on disk, so this figure is
%   a single-epsilon result and the caption should say eps = 0.05.
%
%   Output: figures_out/risk_allocation_ssh.{pdf,png}

clear; clc; close all;

%% ---------------------------------------------------------------- CONFIG
here = fileparts(mfilename('fullpath'));
ROOT = fileparts(here);                        % the repository root
FIGD = fullfile(ROOT,'figures_out');
if ~isfolder(FIGD), mkdir(FIGD); end

THETA   = 400;                                 % headline tracking weight
EPSTAG  = 'e050';                              % = SaveTrajEpsilon 0.05
SCENS   = 1:4;
NBINS   = 10;                                  % slack deciles
MCM     = 1.23348;                             % 10^6 m^3 per kaf
FS      = 14;

UNITS = {'McNary','John Day','The Dalles','Bonneville'};
SHORT = {'MCN','JDA','TDA','BON'};             % report columns only
C = [  0 114 178;                              % CVD-safe, checked
       0 158 115;
     213  94   0;
     106  61 154]/255;
REDL = [192  58  43]/255;                      % eps/n reference

%% ------------------------------------------------------------ LOAD DATA
alloc = cell(numel(SCENS),1);
slack = cell(numel(SCENS),1);
for k = 1:numel(SCENS)
    f = fullfile(ROOT,'results','runs', ...
                 sprintf('BPA_HL2_TH%07d',THETA),'jobs', ...
                 sprintf('s%02d_%s_DDU_k100_SSH.mat',SCENS(k),EPSTAG));
    assert(isfile(f), 'missing run: %s', f);
    S = load(f,'traj');
    alloc{k} = S.traj.risk_allocation;                       % (T x n)
    slack{k} = MCM * min(S.traj.V - S.traj.qLow, [], 2);     % (T x 1), 10^6 m^3
end
n = size(alloc{1},2);
assert(n == numel(UNITS), 'expected %d units, found %d', numel(UNITS), n);

pooled = vertcat(slack{:});
fprintf('%d scenarios, %d hours each, theta = %d, eps = 0.05\n', ...
        numel(SCENS), size(alloc{1},1), THETA);

%% ------------------------------------------------- bin by slack decile
edg  = quantile(pooled, linspace(0,1,NBINS+1));
edg(1) = -inf; edg(end) = inf;                 % keep the tails inside

ctr  = nan(NBINS,1);
M    = nan(NBINS,n);                           % mean share, across scenarios
SD   = nan(NBINS,n);                           % scenario spread
pEven= nan(NBINS,1);                           % % of hours exactly at eps/n

for b = 1:NBINS
    per = nan(numel(SCENS),n);
    binSlack = []; binEven = [];
    for k = 1:numel(SCENS)
        m = slack{k} >= edg(b) & slack{k} < edg(b+1);
        per(k,:) = mean(alloc{k}(m,:), 1, 'omitnan');
        binSlack = [binSlack; slack{k}(m)];                       %#ok<AGROW>
        binEven  = [binEven;  all(abs(alloc{k}(m,:) - 1/n) < 1e-6, 2)]; %#ok<AGROW>
    end
    ctr(b)   = median(binSlack);
    M(b,:)   = mean(per,1);
    SD(b,:)  = std(per,0,1);
    pEven(b) = 100*mean(binEven);
end

% Regime bounds read off the data, not chosen by eye: the dynamic band ends
% where fewer than 10% of hours are exactly even, the even band starts where
% more than 90% of them are.
bDyn  = find(pEven <  10, 1, 'last' );
bEven = find(pEven >= 90, 1, 'first');
if isempty(bDyn),  bDyn  = 0;      end         % no clearly dynamic decile
if isempty(bEven), bEven = NBINS+1; end        % never settles onto eps/n
qDyn  = quantile(pooled, bDyn/NBINS);
qEven = quantile(pooled, (bEven-1)/NBINS);

%% ---------------------------------------------------------------- FIGURE
fig = figure('Units','inches','Position',[1 1 7.16 3.8],'Color','w');
ax  = axes(fig); hold(ax,'on');

xl = [min(ctr)-2, max(ctr)+3];
yl = [0 0.51];

grid(ax,'on');
ax.YGrid = "on";
ax.XGrid = "off";
% plot() rather than yline(): yline reads a trailing char arg as a label, so
% yline(y,'--','Color',...) misparses on some releases
plot(ax, xl, [1 1]/n, '--', 'Color', "k", 'LineWidth', 1.4, ...
     'HandleVisibility','off');

for i = 1:n
    lo = M(:,i) - SD(:,i);  hi = M(:,i) + SD(:,i);
    patch(ax, [ctr; flipud(ctr)], [lo; flipud(hi)], C(i,:), ...
          'FaceAlpha',0.16, 'EdgeColor','none', 'HandleVisibility','off');
end
for i = 1:n
    plot(ax, ctr, M(:,i), '-o', 'Color', C(i,:), 'LineWidth', 1.7, ...
         'MarkerSize', 4.6, 'MarkerFaceColor', C(i,:), ...
         'MarkerEdgeColor', C(i,:), 'DisplayName', UNITS{i});
end

xlim(ax, [min(ctr), max(ctr)]); ylim(ax, yl);
xlabel(ax, 'Volume Buffer to the P10 floor (10^6 m^3)');
ylabel(ax, 'Risk Share $\varepsilon_{i,t}$', ...
       'Interpreter','latex');

% eps/n replaces the 0.25 tick, so the reference needs no annotation
yticks(ax, [0 0.1 0.2 1/n 0.3 0.4 0.5 0.6]);
yticklabels(ax, {'0','0.1','0.2','\epsilon/n','0.3','0.4','0.5','0.6'});


legend(ax, 'Location','southoutside', 'Box','off', 'NumColumns', 4);
ax.FontSize = FS; ax.FontName = 'Times New Roman';
ax.LineWidth = 0.7; ax.XColor = [0.25 0.25 0.25]; ax.YColor = [0.25 0.25 0.25];
ax.GridColor = [0.86 0.86 0.85]; ax.GridAlpha = 1; ax.Box = 'on';

exportgraphics(fig, fullfile(FIGD,'risk_allocation_ssh.pdf'), 'ContentType','vector');
exportgraphics(fig, fullfile(FIGD,'risk_allocation_ssh.png'), 'Resolution', 400);

%% ---------------------------------------------------------------- REPORT
allA  = vertcat(alloc{:});
isEven= all(abs(allA - 1/n) < 1e-6, 2);
nEff  = 1 ./ sum(allA.^2, 2);
[~,who] = max(allA(~isEven,:), [], 2);

fprintf('regime bands: dynamic below %.1f, even above %.1f (10^6 m^3)\n', qDyn, qEven);
fprintf('exactly even in %.1f%% of hours; effective units sharing the budget %.2f even vs %.2f dynamic\n', ...
        100*mean(isEven), mean(nEff(isEven)), mean(nEff(~isEven)));
fprintf('when concentrated, the budget lands on:\n');
for i = 1:n
    fprintf('   %-12s %5.1f%%\n', UNITS{i}, 100*mean(who == i));
end
fprintf('\n slack (10^6 m^3)   %%even   %s\n', sprintf('%8s', SHORT{:}));
for b = 1:NBINS
    fprintf('   %8.1f      %6.1f   %s\n', ctr(b), pEven(b), sprintf('%8.3f', M(b,:)));
end
fprintf('wrote risk_allocation_ssh.{pdf,png} to %s\n', FIGD);

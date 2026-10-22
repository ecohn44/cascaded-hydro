%MAKE_KAPPA_FIGURE  Generation and violation against the DDU risk attitude.
%
%   Run it:  >> make_kappa_figure        (script, not a function)
%
%   Replaces the old Fig. 9, which plotted both quantities against a raw
%   upstream release coefficient on a 0.005-0.03 axis. That axis cannot be
%   reused: the fitted coefficients in the four-project model are 0.087
%   (MCN), 0.189 (JDA), 0.601 (TDA) and 0.185 (BON), one to twenty times
%   larger. The sweep is therefore over kappa, the multiplier on the fitted
%   ramp-dependent slope, with kappa = 1 the fitted model.
%
%   WHAT KAPPA IS. dispatchSetup.m sets gammaModel = cfg.Kappa*gammaData and
%   then rebalances the intercept, sig0 = max(sigDIU - gamma*meanDu,
%   0.2*sigDIU), so the fitted conditional standard deviation AT THE MEAN
%   RAMP is preserved. Kappa is a pure SHAPE parameter -- how much variance
%   the controller attributes to ramping rather than to baseline -- not a
%   conservatism dial. The realised world never moves with it: the shock
%   path is built from gammaData, never gammaModel. That is why both panels
%   can improve at once, and it is worth saying in the text, because a
%   reviewer will ask why more caution is free.
%
%   KAPPA = 0 IS THE DIU CONTROL. dispatchSetup takes the identical arithmetic
%   path at Kappa = 0, so the leftmost point of each panel is the
%   decision-independent policy: generation increase is zero by
%   construction and IVI is the DIU level. Say so in the caption.
%
%   THE GRID STOPS AT 1.5. The sig0 clamp binds when
%       kappa > 0.8*sigDIU/(gammaData*meanDu)
%   which on the WY2024 fit is 16.69 (MCN), 8.28 (JDA), 6.22 (BON) and
%   1.52 (TDA). The Dalles binds first because 52.6% of its variance is
%   already ramp-attributed against 5-13% elsewhere. Past 1.52 its
%   intercept sits on the floor and further kappa RAISES its total sigma
%   instead of redistributing it, so the sweep stops measuring risk
%   attitude. Do not extend the grid without re-deriving that bound.
%
%   Each point is the mean over the four inflow scenarios, and generation
%   is differenced against each scenario's OWN kappa = 0 run, so the
%   comparison is paired. All four scenarios favour DDU at every kappa on
%   both metrics. Bands are deliberately not drawn: between-scenario IVI
%   levels span 868 to 3280 (10^6 m^3), far wider than the effect, so a
%   min-max band would cross the DIU level and imply DDU is sometimes
%   worse, which the paired comparison says it never is.
%
%   NOTE the kappa = 1.00 point sits above its neighbours on IVI, 1650
%   against 1497 at 0.75. That is scenario noise, not structure: the paired
%   per-scenario differences against kappa = 0.75 are -1.1, -14.8, +3.7 and
%   -10.8 percentage points, t = -1.4 on 3 d.o.f., p = 0.26. A clause in
%   the caption stops a reader reading a mechanism into it.
%
%   Source: results/runs/BPA_KAP_TH0000400_E050, 28 jobs,
%   theta = 400, eps = 0.05, SSH, four scenarios, WY2024.
%
%   Output: figures_out/kappa_sweep.{pdf,png}

clear; clc; close all;

%% ---------------------------------------------------------------- CONFIG
here = fileparts(mfilename('fullpath'));
ROOT = fileparts(here);                        % the repository root
FIGD = fullfile(ROOT,'figures_out');
if ~isfolder(FIGD), mkdir(FIGD); end

TAG = 'BPA_KAP_TH0000400_E050';
MCM = 1.23348;                                 % 10^6 m^3 per kaf
FS  = 16;

RED  = [0.8500 0.3250 0.0980];                 % generation, as the old figure
DARK = [0.25   0.25   0.25  ];                 % violation, as the old figure

%% ------------------------------------------------------------ LOAD DATA
f = fullfile(ROOT,'results','runs',TAG,sprintf('summary_%s.csv',TAG));
assert(isfile(f), 'missing summary: %s', f);
T = readtable(f);

kappas = unique(round(T.kappa,2))';
scens  = unique(T.scen)';
nK = numel(kappas); nS = numel(scens);

dGen = nan(nK,nS);     % GWh above this scenario's own kappa = 0 run
ivi  = nan(nK,nS);     % 10^6 m^3
for a = 1:nK
    for b = 1:nS
        r = T(abs(round(T.kappa,2)-kappas(a))<1e-9 & T.scen==scens(b), :);
        c = T(abs(T.kappa)<1e-9                    & T.scen==scens(b), :);
        assert(height(r)==1 && height(c)==1, ...
            'expected one row for kappa %.2f scenario %d', kappas(a), scens(b));
        dGen(a,b) = (r.energy_MWh - c.energy_MWh)/1e3;
        ivi(a,b)  = r.lower_ivi_kaf * MCM;
    end
end
mGen = mean(dGen,2);
mIvi = mean(ivi, 2);

%% ---------------------------------------------------------------- FIGURE
fig = figure('Units','inches','Position',[1 1 8 4],'Color','w');
tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'off');
ax1.GridAlpha = 0.25;
plot(ax1, kappas, mGen, '-o', 'Color',RED, 'LineWidth',2.5, ...
     'MarkerSize',7, 'MarkerFaceColor',RED, 'MarkerEdgeColor',RED);
ylabel(ax1, {'Mean Generation','Increase (GWh)'});
ylim(ax1, [-0.8 16.5]);

ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'off');
ax2.GridAlpha = 0.25;
plot(ax2, kappas, mIvi, '-o', 'Color',DARK, 'LineWidth',2.5, ...
     'MarkerSize',7, 'MarkerFaceColor',DARK, 'MarkerEdgeColor',DARK);
ylabel(ax2, {'Mean IVI','(10^6 m^3)'});
xlabel(ax2, 'Decision-Dependence Scaling (\kappa)');
ylim(ax2, [1280 2200]);

for ax = [ax1 ax2]
    xlim(ax, [0 1.5]);
    xticks(ax, 0:0.25:1.5);
    ax.FontName = 'Times New Roman';  ax.FontSize = FS;
    ax.LineWidth = 0.8;  ax.Box = 'on';  ax.TickDir = 'in';
    ax.XColor = [0.30 0.30 0.30];  ax.YColor = [0.30 0.30 0.30];
    ax.GridColor = [0.85 0.85 0.85];  ax.GridAlpha = 1;
end
xticklabels(ax1, []);

exportgraphics(fig, fullfile(FIGD,'kappa_sweep.pdf'), 'ContentType','vector');
exportgraphics(fig, fullfile(FIGD,'kappa_sweep.png'), 'Resolution', 400);

%% ---------------------------------------------------------------- REPORT
fprintf('%s | %d kappa x %d scenarios\n', TAG, nK, nS);
fprintf('  kappa   dGen (GWh)    IVI (10^6 m^3)   IVI vs DIU   scenarios favouring DDU\n');
for a = 1:nK
    red = 100*(mIvi(1) - mIvi(a))/mIvi(1);
    win = sum(ivi(a,:) < ivi(1,:));
    fprintf('   %4.2f     %+7.2f        %8.0f       %+6.1f%%        %d/%d\n', ...
            kappas(a), mGen(a), mIvi(a), red, win, nS);
end
fprintf('headline: IVI %.0f -> %.0f (10^6 m^3), %.1f%% lower, generation +%.1f GWh\n', ...
        mIvi(1), mIvi(end), 100*(1-mIvi(end)/mIvi(1)), mGen(end));
fprintf('wrote kappa_sweep.{pdf,png} to %s\n', FIGD);

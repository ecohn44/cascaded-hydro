%MAKE_FIG2  Fig. 2 for the TSTE revision -- four SEPARATE panel figures.
%
%   Run it:  >> make_fig2          (script, not a function -- every variable
%   is left in the base workspace so you can poke at the numbers and restyle
%   panels one at a time. The helpers at the bottom are local functions,
%   which a script is allowed to have since R2016b.)
%
%   v3 (2026-09-26): each panel is now its own figure/file, with no in-figure
%   titles -- panel titles/captions belong in the LaTeX \subfloat/\subcaption
%   wrapper instead. Axes carry only clear, capitalized x/y labels (plus, for
%   panel (b), a small in-axes dam-name label, since that identifies which
%   of the 2x2 quadrant is which dam -- a caption on the outer subfigure
%   can't reach inside it).
%     (a) fig2a_hydrograph        Lower Columbia hydrograph, four-project
%                                 mean, p10/mean/p90, Jan-Dec. Unchanged.
%     (b) fig2b_routing           AR mean model on held-out years, observed
%                                 vs. predicted density, one small panel per
%                                 dam (McNary, John Day, The Dalles,
%                                 Bonneville), kept as a single 2x2 image.
%     (c) fig2c_forecast_error    ARCH-X's fitted Normal(0,sigma) overlaid
%                                 directly (one curve per release bin, color
%                                 = release), with the constant DIU sigma
%                                 shown as a reference curve. Skew is
%                                 reported separately (results table), not
%                                 in this panel -- it shows heteroskedasticity
%                                 only, matching how the uncertainty model is
%                                 actually used (Gaussian chance constraint).
%     (d) fig2d_ddu_window        A real 2-week window where DDU's ARCH-X
%                                 band widens with the JDA ramp and catches
%                                 the forecast-error spikes that a constant
%                                 DIU band would flag as violations.
%     (e) fig5_forecast_error_release  Std. dev. of TDA forecast error vs.
%                                 effective upstream release (test years):
%                                 binned empirical points, the constant DIU
%                                 line, and the smooth ARCH-X-fitted DDU
%                                 curve. Same 6-bin scheme as panel (c).
%
%   Reads only the raw USACE exports. Every fitted coefficient (mean model +
%   ARCH-X, all 4 dams) is hard-coded in the CONFIG block below, taken from
%       <desk>/cascaded-hydro/streamflow-models/fit_inflow_fullyear.py
%   (spill-fixed, ridge24 mean specification, train WY2018-22 / test
%   WY2023-25). Change a number there and this script redraws; nothing is
%   refitted here, so MATLAB cannot drift from the Python fit.
%
%   Outflow is TOTAL release, Flow-Gen + Flow-Spill. This script reads the
%   raw USACE exports only -- it does not touch data/inputs.mat or anything
%   else the controller uses, and it changes no simulation input.
%
%   Output: figures_out/fig2a_hydrograph.{pdf,png}
%           figures_out/fig2b_routing.{pdf,png}
%           figures_out/fig2c_forecast_error.{pdf,png}
%           figures_out/fig2d_ddu_window.{pdf,png}
%           figures_out/fig5_forecast_error_release.{pdf,png}

clear; clc;

%% ---------------------------------------------------------------- CONFIG
% Paths resolve from this file, matching runSweep.m's convention, so the
% script works wherever the project tree is checked out. Expected layout:
%     <desk>/matlab_ddu_four_unit/paper_figures/make_fig2.m   <- this file
%     <desk>/cascaded-hydro/simulation-data/<project>/<project>.csv
%     <desk>/Columbia River Data/US ACE Data 2026/<project>-head.csv
here = fileparts(mfilename('fullpath'));
ROOT = fileparts(here);                       % the repository root
% Raw USACE exports ship with the repository, so nothing is read from
% outside the checkout.
SIM  = fullfile(ROOT, 'data', 'raw', 'usace-flow');
ACE  = fullfile(ROOT, 'data', 'raw', 'usace-forebay');
FIGD = fullfile(ROOT, 'figures_out');

assert(isfolder(SIM), 'make_fig2:simMissing', ...
    'Flow data not found at %s -- edit SIM in the CONFIG block.', SIM);
assert(isfolder(ACE), 'make_fig2:aceMissing', ...
    'Forebay data not found at %s -- edit ACE in the CONFIG block.', ACE);
if ~isfolder(FIGD), mkdir(FIGD); end

units  = {'mcn','jda','tda','bon'};
proj   = {'mcnary','johnday','dalles','bonneville'};
pre    = {'MCN','JDA','TDA','BON'};
headA  = [1.00 0.83 0.73 0.85];          % forebay->volume shape, per unit
headB  = [0.00 0.12 0.21 0.08];
K      = 0.0826446;                      % kaf per kcfs-hour
CMS    = 28.3168;                        % m^3/s per kcfs
TRAIN  = 2018:2022;   TEST = 2023:2025;

% --- per-unit mean model (ridge24, own lags 1:24 + upstream release lags)
% and ARCH-X conditional variance, all copied verbatim from
% fit_inflow_fullyear.py's fullyear_ridge24 spec. mcn has no upstream unit
% in this reach, so its "model" is pure one-lag persistence and its ARCH-X
% has no release/ramp term (matches the Python fit's mcn special-case).
MDL = struct('name',{},'upstream',{},'tau',{},'B0',{},'QLAG',{},'QCOEF',{}, ...
             'RLAG',{},'RCOEF',{},'omega',{},'alpha',{},'g_rel',{},'g_ramp',{}, ...
             'v0',{},'rel_min',{},'rel_rng',{},'rmp_min',{},'rmp_rng',{}, ...
             'sigma_DIU',{});

MDL(1) = struct('name','mcn','upstream','','tau',NaN, ...
    'B0',0.0,'QLAG',1,'QCOEF',1.0,'RLAG',[],'RCOEF',[], ...
    'omega',68.39245547352047,'alpha',0.2656021724283783,'g_rel',0.0,'g_ramp',0.0, ...
    'v0',89.85736664054856,'rel_min',0.0,'rel_rng',1.0,'rmp_min',0.0,'rmp_rng',1.0, ...
    'sigma_DIU',9.479312561602162);

MDL(2) = struct('name','jda','upstream','mcn','tau',10, ...
    'B0',0.24225651729057063,'QLAG',1:24,'QCOEF',[ ...
      1.12074676 -0.23143458 -0.00912926 -0.00469043 -0.00495693 -0.02671501 ...
     -0.02952234 -0.02662365  0.06508141  0.00061729  0.01556349  0.01556629 ...
     -0.01286978 -0.01749457 -0.01029470 -0.05393533  0.02530521  0.01839794 ...
      0.02234139  0.00845066  0.01049564  0.03165186  0.13894447 -0.07860834], ...
    'RLAG',[9 10 11],'RCOEF',[0.01957648 0.00695373 0.00548536], ...
    'omega',39.86571805966577,'alpha',0.2252424666315235,'g_rel',51.787117733902754, ...
    'g_ramp',23.634242128892804,'v0',68.3960856481441, ...
    'rel_min',47.9,'rel_rng',469.5,'rmp_min',0.0,'rmp_rng',74.39999999999998, ...
    'sigma_DIU',8.270192600426189);

MDL(3) = struct('name','tda','upstream','jda','tau',2, ...
    'B0',0.5197505746442914,'QLAG',1:24,'QCOEF',[ ...
      0.89714400 -0.15788527  0.02877379 -0.01188400 -0.01417420 -0.03312762 ...
     -0.01455826 -0.00084958  0.05319990  0.00790728  0.01585623 -0.00911670 ...
      0.00331579 -0.01981391 -0.00272361 -0.05251301  0.02447719  0.01244897 ...
      0.01294586  0.01514867  0.00214110  0.02506731  0.13084201 -0.08391526], ...
    'RLAG',[1 2 3],'RCOEF',[0.36615223 -0.21533754 0.00757301], ...
    'omega',35.93363590660807,'alpha',0.1736814331710667,'g_rel',61.12944733824949, ...
    'g_ramp',14.634112398619068,'v0',60.6423311623629, ...
    'rel_min',50.8,'rel_rng',466.3,'rmp_min',0.0,'rmp_rng',101.09999999999997, ...
    'sigma_DIU',7.787318611843418);

MDL(4) = struct('name','bon','upstream','tda','tau',4, ...
    'B0',0.3063089512114914,'QLAG',1:24,'QCOEF',[ ...
      1.23337297 -0.30097188  0.00436717 -0.00683676 -0.01495385 -0.02616897 ...
      0.00246222 -0.01315426  0.03047101  0.01889004  0.00052087  0.00752304 ...
     -0.00382333 -0.00424788 -0.02447322 -0.00501978  0.01339631  0.00644078 ...
      0.01249213  0.00458486  0.01846936  0.02785399  0.09847428 -0.10171162], ...
    'RLAG',[3 4 5],'RCOEF',[0.02864034 -0.00915402 0.00184186], ...
    'omega',10.521644577201169,'alpha',0.2451478751985036,'g_rel',34.47408216689823, ...
    'g_ramp',0.0,'v0',24.914305665223797, ...
    'rel_min',23.4,'rel_rng',468.6,'rmp_min',0.0,'rmp_rng',118.1, ...
    'sigma_DIU',4.991423210390379);

% --- style ---------------------------------------------------------------
DIU_C = [  0 114 178]/255;    % colourblind-safe pair
DDU_C = [213  94   0]/255;
INK   = [ 26  26  26]/255;
FS_A  = 16;                   % panel (a): bigger for paper readability (cap 16-18)
FS_B  = 14;                   % panel (b): a little bigger than the shared FS
FS_C  = 14;                   % panel (c): bigger fonts, per request
FS_D  = 14;                   % panel (d): bigger axis/legend labels
FS_5  = 14;                   % fig. 5: bigger fonts again, per request
FIGW  = 7.16;                 % inches, IEEE double-column figure width

% Panel (c)'s release-color gradient: the 45-85% slice of matplotlib's
% YlOrBr (lightened toward yellow/orange at the dark end vs. the previous
% 55-100% slice, which bottomed out near-black). 32-row lookup table,
% generated once from Python and pasted here verbatim so MATLAB doesn't
% need any colormap toolbox extras -- interp1 does the rest.
DDU_LUT = [
    0.996078 0.666128 0.219223;
    0.996078 0.650258 0.205198;
    0.996078 0.634387 0.191173;
    0.996078 0.618516 0.177147;
    0.994971 0.597478 0.159493;
    0.988328 0.582345 0.151742;
    0.981684 0.567213 0.143991;
    0.972826 0.547036 0.133656;
    0.966182 0.531903 0.125905;
    0.959539 0.516770 0.118155;
    0.950681 0.496594 0.107820;
    0.944037 0.481461 0.100069;
    0.937393 0.466328 0.092318;
    0.928535 0.446151 0.081984;
    0.919093 0.432018 0.074833;
    0.907282 0.418731 0.068189;
    0.891534 0.401015 0.059331;
    0.879723 0.387728 0.052687;
    0.867912 0.374441 0.046044;
    0.856101 0.361153 0.039400;
    0.840354 0.343437 0.030542;
    0.828543 0.330150 0.023899;
    0.816732 0.316863 0.017255;
    0.800984 0.299146 0.008397;
    0.782745 0.289919 0.008520;
    0.763922 0.281061 0.009258;
    0.738824 0.269250 0.010242;
    0.720000 0.260392 0.010980;
    0.701176 0.251534 0.011719;
    0.676078 0.239723 0.012703;
    0.657255 0.230865 0.013441;
    0.638431 0.222007 0.014179;
];

%% ------------------------------------------------------------ LOAD DATA
fb = table(); fl = table();
for i = 1:4
    h = readtable(fullfile(ACE, [proj{i} '-head.csv']), 'VariableNamingRule','preserve');
    t = datetime(h.("Date Time"), 'InputFormat','dd-MMM-yyyy HH:mm');
    v = h.(sprintf('%s.Elev-Forebay.Inst.1Hour.0.CBT-REV [ft]', pre{i}));
    k = table(t, v, 'VariableNames', {'dt', [units{i} '_fb']});
    k = k(~isnat(k.dt) & ~isnan(k.(2)), :);
    fb = mergeOn(fb, k);

    d = readtable(fullfile(SIM, [proj{i} '.csv']), 'VariableNamingRule','preserve');
    t = datetime(d.("Date Time"));
    gen = d.(sprintf('%s.Flow-Gen.Ave.1Hour.1Hour.CBT-REV [kcfs]',   pre{i}));
    spl = d.(sprintf('%s.Flow-Spill.Ave.1Hour.1Hour.CBT-REV [kcfs]', pre{i}));
    inf = d.(sprintf('%s.Flow-In.Inst.~6Hours.0.RFC-FCST [kcfs]',    pre{i}));
    spl(isnan(spl)) = 0;
    k = table(t, gen + spl, inf, 'VariableNames', ...
              {'dt', [units{i} '_out'], [units{i} '_in']});
    k = k(~isnat(k.dt), :);
    fl = mergeOn(fl, k);
end

%% ------------------------------------------------- (a) hydrograph series
inCols = fl{:, cellfun(@(u) [u '_in'], units, 'uni', 0)};
good   = sum(~isnan(inCols), 2) >= 3;
sysq   = mean(inCols(good, :), 2, 'omitnan') * CMS;
sysdt  = fl.dt(good);
dayn   = dateshift(sysdt, 'start', 'day');
[ud, ~, ix] = unique(dayn);
dq  = accumarray(ix, sysq, [], @mean);
yr  = year(ud);  doy = min(day(ud, 'dayofyear'), 365);
keep = ismember(yr, findFullYears(yr, doy));
dq = dq(keep); doy = doy(keep);
p10 = zeros(365,1); p50 = zeros(365,1); p90 = zeros(365,1);
for dd = 1:365
    v = dq(doy == dd);
    p10(dd) = quantile(v, 0.10);  p50(dd) = mean(v);  p90(dd) = quantile(v, 0.90);
end

%% --------------------------------------- reconstruction (all 4 dams' q)
M = innerjoin(fl, fb, 'Keys', 'dt');
M = sortrows(M, 'dt');
ok = all(~isnan(M{:, cellfun(@(u) [u '_fb'], units, 'uni', 0)}), 2);
M  = M(ok, :);

for i = 1:4
    H  = M.([units{i} '_fb']);
    Hn = (H - min(H)) ./ (max(H) - min(H));
    Vn = min(max((Hn - headB(i)) ./ headA(i), 0), 1);
    Vk = min(Vn) + Vn .* (max(Vn) - min(Vn));
    q  = M.([units{i} '_out']) + [NaN; diff(Vk)] / K;
    q  = q - (mean(q, 'omitnan') - mean(M.([units{i} '_out']), 'omitnan'));
    M.([units{i} '_q']) = q;
end
M.sy = year(M.dt) - double(month(M.dt) <= 4);          % season year

%% ------------------------------- per-unit mean model + ARCH-X, all 4 dams
newSeason = [true; diff(M.sy) ~= 0];
RES = struct();
for i = 1:4
    mdl = MDL(i);
    qU  = M.([units{i} '_q']);

    if isempty(mdl.upstream)
        pred = shiftSeason(qU, 1, newSeason);            % mcn: pure persistence
    else
        n  = height(M);
        X  = nan(n, numel(mdl.QLAG) + numel(mdl.RLAG));
        for h = 1:numel(mdl.QLAG)
            X(:,h) = shiftSeason(qU, mdl.QLAG(h), newSeason);
        end
        uup = M.([mdl.upstream '_out']);
        for j = 1:numel(mdl.RLAG)
            X(:, numel(mdl.QLAG)+j) = shiftSeason(uup, mdl.RLAG(j), newSeason);
        end
        % NaN must survive the floor. MATLAB's max IGNORES NaN, so the
        % obvious max(0, ...) turns every blanked season-boundary row into
        % a prediction of 0 and hands back a residual equal to the whole
        % inflow (~6500 m^3/s at TDA). Those rows are finite, so the
        % ~isnan(resid) test below cannot catch them, and they land in the
        % R2, in the ARCH-X e^2 recursion and in any spread computed from E.
        lin  = mdl.B0 + X * [mdl.QCOEF, mdl.RCOEF]';
        pred = max(0, lin);
        pred(isnan(lin)) = NaN;
    end
    resid = qU - pred;

    % held-out (2023-25) predicted-vs-observed sample, for panel (b)
    sel = ismember(M.sy, TEST) & ~isnan(qU) & ~isnan(pred);
    RES(i).yo = qU(sel)*CMS;  RES(i).yp = pred(sel)*CMS;
    RES(i).R2 = 1 - sum((RES(i).yo-RES(i).yp).^2) / sum((RES(i).yo-mean(RES(i).yo)).^2);

    % ARCH-X sigma path (needed in full for TDA's panels (c)/(d); computed
    % for every unit for consistency -- cheap either way)
    if isempty(mdl.upstream)
        rel = nan(height(M),1); ramp = nan(height(M),1);
    else
        uup  = M.([mdl.upstream '_out']);
        rel  = shiftSeason(uup, mdl.tau, newSeason);
        ramp = abs(shiftSeason([NaN; diff(uup)], mdl.tau, newSeason));
    end
    gd = ~isnan(resid) & ~isnan(rel) & ~isnan(ramp);
    E = resid(gd); R = rel(gd); P = ramp(gd); SY = M.sy(gd); DT = M.dt(gd);
    relN = max((R - mdl.rel_min)/mdl.rel_rng, 0);
    rmpN = max((P - mdl.rmp_min)/mdl.rmp_rng, 0);
    rst  = [true; diff(SY) ~= 0];
    e2   = [0; E(1:end-1).^2];
    s2   = mdl.omega + mdl.alpha*e2 + mdl.g_rel*relN + mdl.g_ramp*rmpN;
    s2(rst) = mdl.v0;
    sig  = sqrt(max(s2, 1e-12));

    RES(i).E = E; RES(i).rel = R; RES(i).ramp = P; RES(i).SY = SY;
    RES(i).DT = DT; RES(i).sig = sig;
end
TDA = RES(3);   % panels (c) and (d) stay TDA/JDA-driven, as in the submitted figure

fprintf('R2 by dam: mcn %.3f  jda %.3f  tda %.3f  bon %.3f\n', ...
        RES(1).R2, RES(2).R2, RES(3).R2, RES(4).R2);

%% ============================================================ (a) FIGURE
f = figure('Units','inches','Position',[1 1 FIGW 2.3],'Color','w');
ax = axes(f); hold(ax,'on'); grid(ax,'off');
fill(ax, [1:365, 365:-1:1], [p10; flipud(p90)], DIU_C, ...
     'FaceAlpha',0.22, 'EdgeColor','none', 'DisplayName','10-90%');
plot(ax, 1:365, p50, 'Color', INK, 'LineWidth', 1.8, 'DisplayName','Mean');
xlim(ax,[1 365]); xticks(ax,[1 60 121 182 244 305]);
xticklabels(ax,{'Jan','Mar','May','Jul','Sep','Nov'});
ylabel(ax,'Inflow (m^3/s)');
box(ax, 'on');   
legend(ax,'Location','northwest','Box','off','FontSize',FS_A);
styleAxes(ax, FS_A);

exportgraphics(f, fullfile(FIGD,'fig2a_hydrograph.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(FIGD,'fig2a_hydrograph.png'), 'Resolution', 400);
fprintf('wrote %s\n', fullfile(FIGD,'fig2a_hydrograph.pdf'));

%% ============================================================ (b) FIGURE
%% ============================================================ (b) FIGURE

f = figure('Units','inches','Position',[1 1 FIGW 5.0],'Color','w');

tl = tiledlayout(f, 2, 2, 'TileSpacing','compact', 'Padding','compact');

unitLabel = {'Unit 1: McNary','Unit 2: John Day', ...
             'Unit 3: The Dalles','Unit 4: Bonneville'};

hObs = gobjects(1); hPred = gobjects(1);

for i = 1:4

    ax = nexttile(tl); hold(ax,'on');

    yo = RES(i).yo;
    yp = RES(i).yp;

    % Extend the density curves across the full distribution
    allY = [yo; yp];
    pad  = 0.08*(max(allY) - min(allY));
    lo   = max(0, min(allY) - pad);
    hi   = max(allY) + pad;
    xs   = linspace(lo, hi, 400);

    fo = ksdensity(yo, xs);
    fp = ksdensity(yp, xs);

    fill(ax, [xs fliplr(xs)], [fo zeros(size(fo))], INK, ...
         'FaceAlpha',0.12, 'EdgeColor','none');

    hO = plot(ax, xs, fo, 'Color', INK, 'LineWidth', 2.0);
    hP = plot(ax, xs, fp, '--', 'Color', DDU_C, 'LineWidth', 2.0);

    if i == 1
        hObs = hO;
        hPred = hP;
    end

    text(ax, 0.93, 0.86, unitLabel{i}, 'Units','normalized', ...
         'HorizontalAlignment','right', ...
         'FontSize',FS_B+.5, ...
         'FontWeight','normal');

    xlabel(ax, 'Inflow (m^3/s)');
    ylabel(ax, 'Density');
    ax.YTick = [];

    styleAxes(ax, FS_B);

end

lg = legend([hObs hPred], {'Observed','Predicted'}, 'Box','off', ...
            'Orientation','horizontal', 'FontSize',FS_B+2);
lg.Layout.Tile = 'south';
exportgraphics(f, fullfile(FIGD,'fig2b_routing.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(FIGD,'fig2b_routing.png'), 'Resolution', 400);
fprintf('wrote %s\n', fullfile(FIGD,'fig2b_routing.pdf'));

%% ============================================================ (c) FIGURE
% ARCH-X's fitted Normal(0,sigma) overlaid directly by release bin, with the
% constant DIU sigma as a reference curve -- skew is quantified separately
% (results table), not here, so this panel shows heteroskedasticity only,
% matching how the uncertainty model is actually used (a Gaussian chance-
% constraint reformulation).
isT    = ismember(TDA.SY, TEST);
C_RAMP = (MDL(3).g_ramp*MDL(3).rel_rng)/(MDL(3).g_rel*MDL(3).rmp_rng);
utFull = TDA.rel + C_RAMP*(TDA.ramp - MDL(3).rmp_min);
Ut     = utFull(isT)*CMS/1000;   % effective release, 10^3 m^3/s
St     = TDA.sig(isT)*CMS;       % ARCH-X sigma, m^3/s
sigDIU = MDL(3).sigma_DIU*CMS;

NBIN  = 6;
edges = unique(quantile(Ut, linspace(0,1,NBIN+1)));
edges(1) = edges(1) - 1e-6;
ctr = []; sigBin = [];
for b = 1:numel(edges)-1
    s = Ut >= edges(b) & Ut < edges(b+1);
    if nnz(s) < 200, continue; end
    ctr(end+1)    = mean(Ut(s)); %#ok<AGROW>
    sigBin(end+1) = mean(St(s)); %#ok<AGROW>
end

% ONE sigma colour scale for panels (c) and (f). (c) draws the MODEL sigma
% and (f) the EMPIRICAL sigma, and the empirical range is the wider of the
% two, so the scale is built to span both -- otherwise the same orange would
% mean two different numbers across the two panels, or (f)'s end bars would
% saturate. This probe repeats (f)'s binning and takes only the spread; the
% bootstrap stays in (f).
Ep = TDA.E(isT)*CMS;
edgesP = unique(quantile(Ut, linspace(0,1,12)));  edgesP(1) = edgesP(1) - 1e-6;
sdProbe = [];
for b = 1:numel(edgesP)-1
    sp = Ut >= edgesP(b) & Ut < edgesP(b+1);
    if nnz(sp) < 100, continue; end
    sdProbe(end+1) = std(Ep(sp)); %#ok<AGROW>
end

cLo = min([sigBin sdProbe]);  cHi = max([sigBin sdProbe]);
cSpan = cHi - cLo;
cLo = cLo - 0.15*cSpan;            % floor below the smallest value, so even
                                   % that one carries visible tint
lutX = linspace(cLo, cHi, size(DDU_LUT,1));
colorFor = @(v) interp1(lutX, DDU_LUT, min(max(v,cLo),cHi));

f  = figure('Units','inches','Position',[1 1 FIGW 3.4],'Color','w');
ax = axes(f); hold(ax,'on'); grid(ax,'off');
xs = linspace(-1000, 1000, 400);

for k = 1:numel(ctr)
    plot(ax, xs, normpdf(xs, 0, sigBin(k)), 'Color', colorFor(sigBin(k)), 'LineWidth', 2.3, ...
         'HandleVisibility','off');
end

xlabel(ax, 'Forecast Error (m^3/s)');
ylabel(ax, 'PDF');
xlim(ax, [-1000 1000]);
box(ax, 'on');   
legend(ax, 'Location','northeast', 'Box','off', 'FontSize', FS_C);
styleAxes(ax, FS_C);

colormap(ax, DDU_LUT);
caxis(ax, [cLo cHi]);  %#ok<CAXIS> -- caxis kept (not clim) for older-MATLAB compatibility
cb = colorbar(ax);
cb.Label.String = 'DDU Standard Deviation (m^3/s)';
cb.Label.FontSize = FS_C;
cb.FontSize = FS_C;
cb.Box = 'off';

exportgraphics(f, fullfile(FIGD,'fig2c_forecast_error.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(FIGD,'fig2c_forecast_error.png'), 'Resolution', 400);
fprintf('wrote %s\n', fullfile(FIGD,'fig2c_forecast_error.pdf'));

%% ============================================================ (d) FIGURE
f  = figure('Units','inches','Position',[1 1 FIGW 2.8],'Color','w');
ax = axes(f); hold(ax,'on'); grid(ax,'off');
w0 = datetime(2024,7,10); w1 = datetime(2024,7,24);
win = TDA.DT >= w0 & TDA.DT <= w1;
tw = TDA.DT(win); Ew = TDA.E(win)*CMS; sw = TDA.sig(win)*CMS;
sD = MDL(3).sigma_DIU*CMS; zc = 1.959964;

try
    fill(ax, [tw; flipud(tw)], [zc*sD*ones(size(tw)); flipud(-zc*sD*ones(size(tw)))], ...
         DIU_C, 'FaceAlpha',0.18, 'EdgeColor','none', 'DisplayName','DIU');
catch
    plot(ax, tw,  zc*sD*ones(size(tw)), ':', 'Color', DIU_C, 'LineWidth',1.5,'DisplayName','DIU');
    plot(ax, tw, -zc*sD*ones(size(tw)), ':', 'Color', DIU_C, 'LineWidth',1.5, 'HandleVisibility','off');
end
plot(ax, tw,  zc*sw, 'Color', DDU_C, 'LineWidth', 1.5, 'HandleVisibility','off');
plot(ax, tw, -zc*sw, 'Color', DDU_C, 'LineWidth', 1.5, 'DisplayName','DDU');
plot(ax, tw, Ew, 'Color', INK, 'LineWidth', 1.2, 'DisplayName','Realized Forecast Error');

yline(ax, 0, 'Color', [0.85 0.85 0.85], 'HandleVisibility','off');
ylabel(ax, 'Forecast Error (m^3/s)');
xlim(ax, [min(tw), max(tw)]);
box(ax, 'on');   
lg = legend(ax, 'Location','southoutside', 'Box','off', 'Orientation','horizontal', ...
            'FontSize', FS_D);
styleAxes(ax, FS_D);

exportgraphics(f, fullfile(FIGD,'fig2d_ddu_window.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(FIGD,'fig2d_ddu_window.png'), 'Resolution', 400);
fprintf('wrote %s\n', fullfile(FIGD,'fig2d_ddu_window.pdf'));

%% ============================================================ (e) FIGURE 5
% Std. dev. of TDA forecast error vs. effective upstream release, test
% years -- same 6-bin scheme as panel (c) for the "Observed" points, the
% constant DIU sigma, and a smooth DDU curve.
%
% The DDU curve is a straight-line least-squares fit (ordinary linear
% regression) THROUGH the "Observed" points themselves -- not an
% interpolation that passes through every dot, and not a separately-derived
% model series (both the closed-form omega/alpha/g_rel formula and a fit to
% the raw ARCH-X sigma path undershot the empirical spread badly, especially
% at the high-release end).
%% ============================================================ (e) FIGURE

% Use more bins for Fig. 5 without changing Fig. 2(c)
NBIN5  = 12;
edges5 = unique(quantile(Ut, linspace(0,1,NBIN5+1)));
edges5(1) = edges5(1) - 1e-6;

ctr5 = [];
sd5  = [];

for b = 1:numel(edges5)-1
    s = Ut >= edges5(b) & Ut < edges5(b+1);

    if nnz(s) < 100
        continue;
    end

    ctr5(end+1) = mean(Ut(s));  %#ok<AGROW>
    sd5(end+1)  = mean(St(s));  %#ok<AGROW>
end

% Trend through the fitted DDU standard deviations
pFit  = polyfit(ctr5, sd5, 1);
xfit  = linspace(min(ctr5), max(ctr5), 200);
dduSd = polyval(pFit, xfit);

f  = figure('Units','inches','Position',[1 1 FIGW 3.4],'Color','w');
ax = axes(f); hold(ax,'on'); grid(ax,'off');

scatter(ax, ctr5, sd5, 55, ...
        'MarkerFaceColor',INK, ...
        'MarkerEdgeColor','none', ...
        'DisplayName','Observed');

yline(ax, sigDIU, '--', ...
      'Color',DIU_C, ...
      'LineWidth',2.2, ...
      'DisplayName','DIU');

plot(ax, xfit, dduSd, ...
     'Color',DDU_C, ...
     'LineWidth',2.2, ...
     'DisplayName','DDU');

xlabel(ax, 'Upstream Release (10^3 m^3/s)');
ylabel(ax, 'Standard Deviation (m^3/s)');

xlim(ax, [min(xfit)-0.2, max(xfit)+0.2]);
yMax = max([sd5(:); sigDIU; dduSd(:)])*1.15;
ylim(ax, [160 300]);
box(ax,'on');
legend(ax,'Location','northwest','Box','off','FontSize',FS_5);
styleAxes(ax,FS_5);
exportgraphics(f, fullfile(FIGD,'fig2e_forecast_error_release.pdf'), 'ContentType','vector');
exportgraphics(f, fullfile(FIGD,'fig2e_forecast_error_release.png'), 'Resolution', 400);
fprintf('wrote %s\n', fullfile(FIGD,'fig2e_forecast_error_release.pdf'));

%% =========================================================== (f) FIGURE
% Std. dev. of TDA forecast error vs. effective upstream release, test
% years -- the EMPIRICAL spread, binned, as bars with a bootstrap interval.
%
% This is the companion to (e) and answers a different question. Panel (e)
% bins St = TDA.sig, the FITTED ARCH-X sigma, and fits a line through those
% points; both of its series are the model. This panel bins the standard
% deviation of the RESIDUALS themselves, which is the quantity the model is
% supposed to reproduce, so the two can be read against each other.
%
% Whiskers are a moving-block bootstrap 95% CI, 2000 resamples of 24 h
% blocks. Blocks rather than i.i.d. because hourly residuals are
% autocorrelated and an i.i.d. bootstrap makes the intervals look tighter
% than they are.
%
% It reuses TDA/Ut/isT from panel (c) above -- no second reconstruction, no
% second copy of the mean model.

rng(20260930);

Ef     = TDA.E(isT)*CMS;                 % signed forecast error, m^3/s
NBINF  = 11;
NBOOTF = 2000;
BLKF   = 24;
edgesF = unique(quantile(Ut, linspace(0,1,NBINF+1)));
edgesF(1) = edgesF(1) - 1e-6;

ctrF = []; sdF = []; loF = []; hiF = []; nF = [];
for b = 1:numel(edgesF)-1
    sel = Ut >= edgesF(b) & Ut < edgesF(b+1);
    if nnz(sel) < 100, continue; end
    xb  = Ef(sel);
    nblk = floor(numel(xb)/BLKF);
    Bk   = reshape(xb(1:nblk*BLKF), BLKF, nblk);
    bsv  = zeros(NBOOTF,1);
    for kk = 1:NBOOTF
        bsv(kk) = std(reshape(Bk(:, randi(nblk, nblk, 1)), [], 1));
    end
    ctrF(end+1) = mean(Ut(sel));          %#ok<AGROW>
    nF(end+1)   = nnz(sel);               %#ok<AGROW>
    sdF(end+1)  = std(xb);                %#ok<AGROW>
    loF(end+1)  = quantile(bsv, 0.025);   %#ok<AGROW>
    hiF(end+1)  = quantile(bsv, 0.975);   %#ok<AGROW>
end
nbF = numel(ctrF);

BARF = [0.725 0.827 0.902];
BARE = [0.357 0.486 0.569];
BLK0 = [0 0 0];

f  = figure('Units','inches','Position',[1 1 FIGW 3.4],'Color','w');
ax = axes(f); hold(ax,'on'); grid(ax,'off');

% bars carry the SAME orange ramp as panel (c), on the shared scale set there
hb = bar(ax, 1:nbF, sdF, 0.68, 'FaceColor','flat', 'EdgeColor',BARE, ...
         'LineWidth',0.7);
hb.CData = cell2mat(arrayfun(@(v) colorFor(v), sdF(:), 'UniformOutput',false));
for b = 1:nbF
    plot(ax,[b b],[loF(b) hiF(b)],'-','Color',BLK0,'LineWidth',0.9);
    plot(ax,b+[-1 1]*0.12,[loF(b) loF(b)],'-','Color',BLK0,'LineWidth',0.9);
    plot(ax,b+[-1 1]*0.12,[hiF(b) hiF(b)],'-','Color',BLK0,'LineWidth',0.9);
end

xlim(ax,[0.4 nbF+0.6]);
ylim(ax,[0 max(hiF)*1.08]);
xticks(ax,1:nbF);
xticklabels(ax, compose('%.1f', ctrF));
xlabel(ax, 'Upstream Release (10^3 m^3/s)');
ylabel(ax, 'DDU Standard Deviation (m^3/s)');
box(ax,'on');
styleAxes(ax,FS_5);
ax.XColor = BLK0; ax.YColor = BLK0;
exportgraphics(f, fullfile(FIGD,'fig2f_sigma_release.pdf'),'ContentType','vector');
exportgraphics(f, fullfile(FIGD,'fig2f_sigma_release.png'),'Resolution',400);
fprintf('wrote %s\n', fullfile(FIGD,'fig2f_sigma_release.pdf'));

fprintf('\n U_ctr      n     sigma      boot 95%% CI\n');
for b = 1:nbF
    fprintf('%6.1f %6d %8.1f   [%6.1f, %6.1f]\n', ...
            ctrF(b), nF(b), sdF(b), loF(b), hiF(b));
end

% ---- contamination check. sigma_DIU for TDA is 220.5 m^3/s, so a residual
% past ~2000 is an order of magnitude outside the model and is far more
% likely a reconstruction artefact than hydrology. Print the worst offenders
% with their timestamps so they can be traced in the raw exports.
DTt  = TDA.DT(isT);
[~, ord] = sort(abs(Ef), 'descend');
fprintf('\nresiduals: n %d   mean %+.2f   sd %.1f   n(|e|>2000) %d\n', ...
        numel(Ef), mean(Ef), std(Ef), sum(abs(Ef) > 2000));
fprintf('largest 15 |residual| and when they occur:\n');
for k = 1:min(15, numel(ord))
    fprintf('   %s   resid %+9.1f   rel %6.2f\n', ...
            datestr(DTt(ord(k)), 'yyyy-mm-dd HH:MM'), Ef(ord(k)), Ut(ord(k)));
end

%% ------------------------------------------------------------- helpers
function T = mergeOn(T, k)
if isempty(T), T = k; else, T = innerjoin(T, k, 'Keys','dt'); end
end

function y = shiftSeason(x, h, newSeason)
%SHIFTSEASON lag x by h, blanking across season boundaries
y = [nan(h,1); x(1:end-h)];
bad = false(numel(x),1);
idx = find(newSeason);
for s = idx', bad(s:min(s+h-1, numel(x))) = true; end
y(bad) = NaN;
end

function ys = findFullYears(yr, doy)
u = unique(yr); ys = [];
for k = 1:numel(u)
    if numel(unique(doy(yr==u(k)))) >= 350, ys(end+1) = u(k); end %#ok<AGROW>
end
end

function styleAxes(ax, FS)
ax.FontSize = FS-1; ax.FontName = 'Times New Roman';
ax.LineWidth = 0.6; ax.XColor = [0.33 0.33 0.33]; ax.YColor = [0.33 0.33 0.33];
ax.GridColor = [0.86 0.86 0.85]; ax.GridAlpha = 1; ax.Box = 'on';
end

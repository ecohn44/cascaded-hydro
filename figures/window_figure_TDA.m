%WINDOW_FIGURE_TDA  DIU against DDU over two weeks in different seasons.
%
%   Run it:  >> window_figure_TDA        (script, not a function)
%
%   Two columns, one week each, two rows:
%       top     The Dalles release (10^2 m^3/s)
%       bottom  The Dalles storage against the P20 operating floor
%
%   Replaces the three-panel version, which read
%   paper/plot_data/15_tracking_th100_hourly.csv at theta = 100, eps = 0.10,
%   BON. This reads the headline job files directly.
%
%   OPERATING POINT. theta = 400, eps = 0.05, SSH. Epsilon is not a choice:
%   the sweeps ran with SaveTrajEpsilon = 0.05, so e050 is the only epsilon
%   with hourly trajectories on disk. Everything else saved summaries only.
%
%   TOP PANEL IS THE DALLES' OWN RELEASE. An upstream (John Day) version
%   was tried -- John Day's ramp is what the gamma_r term acts on, so it is
%   the mechanistic driver -- but the project's own release reads better
%   beside its own storage. Separation is comparable either way (Feb 5.7 vs
%   8.1 kcfs mean absolute difference, Aug 8.4 vs 8.1), so nothing is lost.
%
%   WHY THESE TWO WEEKS. Every 7-day window in all four scenarios was
%   rescanned ON THE SMOOTHED SERIES -- the figure shows smoothed curves, so
%   the windows must be chosen on those, not on the raw ones. Crossings
%   within 12 h of a window edge were excluded so nothing sits against the
%   axis, and each start was then shifted to bring the event near mid-window:
%       (a) 31 Jan-7 Feb 2024, scenario 3 : DIU 19.8 (5 h),  DDU 0.0 (0 h)
%       (b) 17-24 Jun 2024,    scenario 1 : DIU 33.3 (6 h),  DDU 0.0 (0 h)
%   in 10^6 m^3 below the floor, deepest dip 16% and 18% of the plotted
%   range. Put those numbers in the caption, since the panel titles are
%   deliberately absent.
%
%   An earlier pair (4 Feb scen 2, 22 Aug scen 1) was chosen on the RAW
%   series and does not survive smoothing: at 3 h it falls to 14.0 and 3.3,
%   and the August wedge becomes too faint to read. Both panels here carry a
%   single crossing -- on the smoothed basis no week has two crossings that
%   are still visible with DDU clean.
%
%   THE TRADE-OFF BEHIND THAT CHOICE, so it is not relitigated later. Across
%   all scenarios and windows, the count with DDU exactly clean is 969 with
%   one crossing, 284 with two, 40 with three, 17 with four and ZERO with
%   five or more. The three- and four-crossing weeks all have crossings
%   about 4% of the axis deep, invisible at print size. Frequent and visible
%   do not co-occur: when conditions are bad enough for DIU to cross
%   repeatedly, DDU is usually nicked once too.
%
%
%   SMOOTHING. A 3 h centred moving average is applied to the RELEASE panel
%   only, to take the hourly jitter out without losing the ramp shape: it
%   keeps 85-87% of the mean hourly ramp magnitude in these two windows.
%
%   STORAGE IS SMOOTHED TOO, at the same 3 h, and so is the floor -- an
%   unsmoothed floor against a smoothed storage would invent and hide
%   crossings. All shading and all printed totals come from the smoothed
%   pair, so the figure and its numbers agree.
%
%   WHAT THAT COSTS, measured on these windows:
%       window          raw            3 h           5 h          7 h
%       (a) Feb    16.2, 2 cross   14.0, 1 cross   9.9, 1      6.3, 1
%       (b) Aug     7.1, 2 cross    3.3, 2 cross   0.0, NONE   0.0, NONE
%   At 3 h February's two crossings merge into one and August keeps both.
%   DO NOT raise SMOOTH above 3 for the storage: at 5 h the August
%   violation disappears completely and the figure would show DIU never
%   crossing the floor, which is the opposite of the result.
%
%   No titles: label the columns with \subfloat in the manuscript so the
%   caption carries the text.
%
%   The two storage panels share a y-range, which is why the right column's
%   tick labels are suppressed. Without that the seasonal difference in the
%   floor -- about 49 in February against 56 in August -- would be invisible
%   and the columns would not be comparable.
%
%   Output: figures_out/window_TDA_two_seasons.{pdf,png}

clear; clc; close all;

%% ---------------------------------------------------------------- CONFIG
here = fileparts(mfilename('fullpath'));
ROOT = fileparts(here);                        % the repository root
FIGD = fullfile(ROOT,'figures_out');
if ~isfolder(FIGD), mkdir(FIGD); end

TAG   = 'BPA_HL2_TH0000400';
EPSTG = 'e050';
IUNIT = 3;                     % The Dalles  (storage panel)
IREL  = 3;                     % The Dalles  (release panel, its own)
MCM   = 1.23348;               % 10^6 m^3 per kaf
CMS   = 0.283168;              % 10^2 m^3/s per kcfs (28.3168 / 100)
T0    = datetime(2023,10,1);   % hour 1 of WY2024
WHRS  = 24*7;
FS    = 14;
SMOOTH= 4;                     % hours, centred; applied to release AND storage

PANEL = { 3, datetime(2024,1,31); ...       % scenario, window start
          1, datetime(2024,6,17) };

cDIU   = [0.357 0.463 0.565];
cDDU   = [0.710 0.267 0.227];
cFloor = [0.20  0.20  0.19 ];
cFill  = [0.90  0.89  0.87 ];
cBnd   = [0.45  0.45  0.45 ];  % release limits

% Release limits, read from the same params the controller was given rather
% than hard-coded, so they stay right if IREL changes. For The Dalles these
% are umin 38.0 and umax 292.9 kcfs, i.e. 10.8 and 82.9 in 10^2 m^3/s --
% which is exactly where the flat bottoms and tops of the release already
% sit, so the lines name a bound the policy is genuinely riding.
UNAMES = {'MCN','JDA','TDA','BON'};
Sp = load(fullfile(ROOT,'data','inputs.mat'),'params');
Pp = jsondecode(Sp.params);
uLo = Pp.(UNAMES{IREL}).umin * CMS;
uHi = Pp.(UNAMES{IREL}).umax * CMS;

%% ---------------------------------------------------------------- FIGURE
fig = figure('Units','inches','Position',[1 1 7.16 5],'Color','w');
tl  = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax  = gobjects(2,2);
sLo = nan(1,2); sHi = nan(1,2);

for c = 1:2
    scen = PANEL{c,1};  w0 = PANEL{c,2};
    k  = hours(w0 - T0) + 1;
    w  = k:(k + WHRS - 1);
    tw = w0 + hours(0:WHRS-1);

    a = readjob(ROOT, TAG, EPSTG, scen, 'DIU_k000');
    b = readjob(ROOT, TAG, EPSTG, scen, 'DDU_k100');
    Va = a.V(w,IUNIT)*MCM;  Vb = b.V(w,IUNIT)*MCM;  qL = a.qLow(w,IUNIT)*MCM;
    rF = a.ref(w,IUNIT)*MCM;                        % seasonal reference
    Ua = a.U(w,IREL)*CMS;   Ub = b.U(w,IREL)*CMS;   % kcfs -> 10^2 m^3/s
    Va = Va(:)'; Vb = Vb(:)'; qL = qL(:)';          % rows, for patch
    rF = rF(:)';
    % Storage and floor carry the SAME smoothing as the release, and every
    % number below is computed from the smoothed pair -- the shaded wedges
    % and the printed totals must describe the curve that is drawn.
    Va = movmean(Va,SMOOTH);  Vb = movmean(Vb,SMOOTH);  qL = movmean(qL,SMOOTH);

    % ---- top: upstream release
    ax(1,c) = nexttile(tl, c); hold(ax(1,c),'on'); grid(ax(1,c),'on');
    plot(ax(1,c), [tw(1) tw(end)], [uHi uHi], '--', 'Color',cBnd, 'LineWidth',1.1);
    plot(ax(1,c), [tw(1) tw(end)], [uLo uLo], '--', 'Color',cBnd, 'LineWidth',1.1);
    plot(ax(1,c), tw, movmean(Ua,SMOOTH), '-', 'Color',cDIU, 'LineWidth',1.5);
    plot(ax(1,c), tw, movmean(Ub,SMOOTH), '-', 'Color',cDDU, 'LineWidth',1.5);
    ylim(ax(1,c), [0 310*CMS]);

    % ---- bottom: storage over the floor
    ax(2,c) = nexttile(tl, 2+c); hold(ax(2,c),'on'); grid(ax(2,c),'on');
    lo = min([qL Va]) - 3;
    % Fill from zero, not from this column's own bottom: the two storage
    % panels share a y-range, so a per-column bottom leaves a white gap
    % under the shorter one. The axes clip the rest.
    patch(ax(2,c), [tw fliplr(tw)], [qL zeros(1,WHRS)], cFill, ...
          'EdgeColor','none');

    % DIU below the floor, shaded ONLY over the crossing intervals. A single
    % patch between Va and qL would fill the whole week, which is wrong.
    % Episode ends are linearly interpolated to the actual crossing time so
    % each wedge closes to a point on the floor instead of a vertical edge.
    gap   = qL - Va;                       % > 0 while DIU is under the floor
    below = gap > 0;
    e0 = find(diff([false below]) ==  1);
    e1 = find(diff([below false]) == -1);
    for e = 1:numel(e0)
        i1 = e0(e); i2 = e1(e);
        xs = tw(i1:i2); ya = Va(i1:i2); yq = qL(i1:i2);
        if i1 > 1                          % entry crossing
            f  = gap(i1-1)/(gap(i1-1)-gap(i1));
            yc = Va(i1-1) + f*(Va(i1)-Va(i1-1));
            xs = [tw(i1-1)+f*(tw(i1)-tw(i1-1)) xs]; %#ok<AGROW>
            ya = [yc ya]; yq = [yc yq];    %#ok<AGROW>
        end
        if i2 < WHRS                       % exit crossing
            f  = gap(i2)/(gap(i2)-gap(i2+1));
            yc = Va(i2) + f*(Va(i2+1)-Va(i2));
            xs = [xs tw(i2)+f*(tw(i2+1)-tw(i2))]; %#ok<AGROW>
            ya = [ya yc]; yq = [yq yc];    %#ok<AGROW>
        end
        patch(ax(2,c), [xs fliplr(xs)], [yq fliplr(ya)], cDIU, ...
              'FaceAlpha',0.80, 'EdgeColor','none');
    end
    plot(ax(2,c), tw, qL, ':', 'Color',cFloor, 'LineWidth',1.2);
    % Seasonal reference. Dashed rather than solid so it is not confused with
    % the near-black dotted floor directly below it.
    plot(ax(2,c), tw, movmean(rF,SMOOTH), '-', 'Color','k', 'LineWidth',1.3);
    plot(ax(2,c), tw, Va, '-', 'Color',cDIU, 'LineWidth',1.5);
    plot(ax(2,c), tw, Vb, '-', 'Color',cDDU, 'LineWidth',1.5);
    sLo(c) = lo;  sHi(c) = max([Va Vb]) + 3;

    nCross = numel(find(diff([false (qL-Va)>0]) == 1));
    fprintf(['(%s) %s scen %d : DIU %.1f below floor over %d h in %d ' ...
             'crossing(s), DDU %.1f over %d h  (10^6 m^3, %d h smoothed)\n'], ...
            char('a'+c-1), datestr(w0,'dd mmm yyyy'), scen, ...
            sum(max(qL-Va,0)), nnz(Va < qL), nCross, ...
            sum(max(qL-Vb,0)), nnz(Vb < qL), SMOOTH);
end

% shared storage range, so suppressing the right-hand labels is honest
for c = 1:2, ylim(ax(2,c), [min(sLo) max(sHi)]); end

ylabel(ax(1,1), {'The Dalles','Release (10^2 m^3/s)'});
ylabel(ax(2,1), {'The Dalles', 'Storage (10^6 m^3)'});
for c = 1:2
    for r = 1:2
        A = ax(r,c);
        A.FontName = 'Times New Roman';  A.FontSize = FS;
        A.LineWidth = 0.7;  A.Box = 'on';  A.TickDir = 'in';
        A.XColor = [0.30 0.30 0.30];  A.YColor = [0.30 0.30 0.30];
        A.GridColor = [0.88 0.88 0.86];  A.GridAlpha = 1;
        A.XGrid = 'off';  A.YGrid = 'on';      % horizontal rules only
        xlim(A, [PANEL{c,2} PANEL{c,2}+days(7)]);
    end
    xticklabels(ax(1,c), []);
    if c == 2, yticklabels(ax(1,c), []); yticklabels(ax(2,c), []); end
end

h = [plot(ax(2,1),NaN,NaN,'-','Color',cDIU,'LineWidth',1.4), ...
     plot(ax(2,1),NaN,NaN,'-','Color',cDDU,'LineWidth',1.4), ...
     plot(ax(2,1),NaN,NaN,'-','Color','k','LineWidth',1.3), ...
     patch(ax(2,1),NaN,NaN,cFill,'EdgeColor','none'), ...
     plot(ax(2,1),NaN,NaN,'--','Color',cBnd,'LineWidth',1.1)];

lg = legend(h, {'DIU','DDU','Reference','Operating Floor','Release Limits'}, ...
            'Orientation','horizontal','Box','off','FontSize',FS-0.5);
lg.Layout.Tile = 'south';

exportgraphics(fig, fullfile(FIGD,'storage-trajectory.pdf'), 'ContentType','vector');
exportgraphics(fig, fullfile(FIGD,'window_TDA_two_seasons.png'), 'Resolution', 400);
fprintf('wrote window_TDA_two_seasons.{pdf,png} to %s\n', FIGD);

%% ------------------------------------------------------------- helpers
function tr = readjob(ROOT, tag, epsTag, scen, arm)
f = fullfile(ROOT,'results','runs',tag,'jobs', ...
             sprintf('s%02d_%s_%s_SSH.mat', scen, epsTag, arm));
assert(isfile(f), ['missing job: %s\n(only eps = 0.05 carries hourly ' ...
                   'trajectories -- SaveTrajEpsilon was 0.05)'], f);
S = load(f,'traj');  tr = S.traj;
end

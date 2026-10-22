%MAKE_SSH_VS_BON  Bonferroni against supporting hyperplanes, same state.
%
%   Run it:  >> make_ssh_vs_bon        (script, not a function)
%
%   Replaces the old Fig. 11, which plotted mean total energy against a
%   drought-intensity multiplier for SSH and BON in two flow regimes. That
%   figure compared two SEPARATE trajectories, so the gap it showed mixed
%   the method difference with everything the two runs did differently
%   afterwards -- and at the headline operating point the annual energy
%   difference is 0.3 GWh out of 26 600, which is 0.001% and invisible.
%   Annual aggregates are the wrong instrument here.
%
%   WHAT THIS FIGURE USES INSTEAD. realtimeSimulation.m already runs a same-state
%   probe (lines 791-820): inside every SSH run it re-solves the Bonferroni
%   model each hour from the IDENTICAL storage state, identical conditional
%   covariance and identical price path, and records
%       ProbeBonFeasible(t)   Bonferroni admits a plan with no storage slack
%       SSHStepFeasible(t)    SSH attains 1 - eps at every look-ahead hour
%   BOTH CURVES HERE COME FROM THE SSH RUN. They are two surrogates
%   evaluated on one trajectory, not two trajectories, which is the only
%   way to attribute the difference to the method. Say this in the caption
%   or a reviewer will assume two runs and discount the whole figure.
%
%   AXES. The x axis is the joint reliability level 1 - eps in percent, 75
%   to 99, to match the monthly violation heatmaps; the underlying sweep is
%   eps in {0.25, 0.20, 0.15, 0.10, 0.05, 0.01}, which lands exactly on
%   {75, 80, 85, 90, 95, 99}. Note the grid is NOT uniform -- the step from
%   95 to 99 is four points, the rest are five -- so the points are plotted
%   at their true positions and the ticks are set explicitly. Do not
%   substitute a categorical axis; the spacing carries information, since
%   the curves bend between 95 and 99. The y axis is raw hours per water
%   year out of 8760, not a percentage, so the band can be read directly as
%   a count of decisions, and it counts hours the target is MET rather than
%   missed, so SSH sits above Bonferroni and higher reads as better. The
%   probe records misses, so the plotted series are complements; that is a
%   reflection, not a rescaling, and the band is unchanged by it.
%
%   WHY THE GAP BETWEEN THE CURVES IS EXACTLY THE RECOVERED HOURS. The SSH
%   feasible set contains the Bonferroni set, so SSH cannot be infeasible
%   at an hour where Bonferroni is feasible. Measured: 0 counterexamples in
%   420 480 probed hours. Containment therefore makes SSH-infeasible a
%   strict SUBSET of BON-infeasible, and
%       BONinfeasible - SSHinfeasible == probe_ssh_only_feasible_steps
%   identically, not approximately: 780.8 - 633.8 = 147.0 at eps = 0.05.
%   The script asserts this at every point. That identity is why one panel
%   with two lines says everything the earlier two-panel version said.
%
%   THE SECOND CLAIM, WHICH IS A ZERO AND SO IS NOT PLOTTED. On hours where
%   both are feasible SSH's look-ahead objective is never worse: 0 hours
%   out of 420 480. Without it a reader can dismiss the gap as SSH merely
%   relaxing the constraint. With it, the gap is a strict gain. It belongs
%   in the caption. The mean objective improvement on the commonly feasible
%   hours (0.79% at 95%, rising to 1.19% at 75%) is printed below for the
%   paragraph -- it is a paired difference and cannot be drawn as two
%   lines, so it is a sentence, not a panel.
%
%   WHY BOTH CURVES RISE WITH RELIABILITY, which looks backwards and will
%   be asked about: a tighter target is met MORE often, not less. At a
%   fixed state a tighter target is of course harder. But a policy held to
%   a tighter target keeps more storage headroom, so the states it actually
%   visits are slacker, and the second effect wins. This is a property of
%   where the trajectory goes, not of either surrogate, and it is common to
%   both curves, so it does not touch the comparison. What DOES matter is
%   that the gap WIDENS as reliability is relaxed: Bonferroni falls from
%   8087 h met at 99% to 7869 h at 75% while SSH holds near 8110 h from 90%
%   down, because SSH is approaching the hydrological floor of hours at
%   which no allocation of eps whatsoever suffices. Bonferroni keeps paying
%   for its fixed split; SSH stops.
%
%   WHY THE MECHANISM IS ALLOCATION, NOT THE UNION BOUND, and say this in
%   the text or a reviewer will assume otherwise. Under independence the
%   Bonferroni-to-exact gap in the scalar multiplier is essentially
%   CONSTANT in n -- 0.0074 at n = 4 against 0.0079 at n = 16. Making the
%   cascade longer does not make Bonferroni worse. What SSH exploits is
%   that the four projects are heterogeneous, so the efficient split of eps
%   across them is nowhere near uniform; that is Fig. 9. This figure is the
%   consequence of Fig. 9 and the two should be read together.
%
%   DDU ONLY. The same curves under DIU sit lower and tell the identical
%   story (SSH recovers 13.1% of Bonferroni's infeasible hours at 99%
%   rising to 22.0% at 75%, against 14.6% and 26.2% here), but DIU is not
%   the paper's operating point and a second pair of lines invites the
%   reader to compare DIU with DDU, which is a different question answered
%   by Fig. 5. Set POLIC below if a reviewer asks whether the advantage
%   survives a decision-independent sigma. It does.
%
%   THE SCENARIO SPREAD IS NOT DRAWN as a band. The four scenarios agree on
%   the ordering and the slope at every point and the spread is roughly
%   plus or minus 15% of the recovered count; a band would obscure a gap
%   that is itself only 100-230 hours. Per-scenario values print below.
%
%   THE COST, for the paragraph rather than the figure: 186 ms per rolling
%   hour against 32 ms for Bonferroni, a factor of 5.8. Both are four to
%   five orders of magnitude inside a 3600 s decision cycle, so the ratio
%   is the honest number to quote and the headroom is what makes it
%   irrelevant. Note the ratio GROWS with cascade size (5.5x at n = 4,
%   12.7x at n = 8, 24.7x at n = 12 in run_scalability), so quote the ratio
%   at n = 4 here and let the scalability table carry the trend.
%
%   Source: results/runs/BPA_HL2_TH0000400, theta = 400,
%   WY2024, 6 epsilons x {DIU, DDU} x {BON, SSH} x 4 scenarios = 96 runs.
%
%   Output: figures_out/ssh_vs_bon.{pdf,png}

clear; clc; close all;

%% ---------------------------------------------------------------- CONFIG
here = fileparts(mfilename('fullpath'));
ROOT = fileparts(here);                        % the repository root
FIGD = fullfile(ROOT,'figures_out');
if ~isfolder(FIGD), mkdir(FIGD); end

TAG   = 'BPA_HL2_TH0000400';
POLIC = 'DDU';                                 % the paper's operating point
FS    = 9.5;

cBON  = [0      0.4470 0.7410];                % blue,  as the original Fig. 11
cSSH  = [0.3600 0.5250 0.1450];                % green, darkened from the
                                               % original [0.4660 0.6740
                                               % 0.1880], which sat too light
                                               % against the fill band
cFill = [0.92   0.91   0.88  ];                % the recovered-hours band

%% ------------------------------------------------------------ LOAD DATA
f = fullfile(ROOT,'results','runs',TAG,sprintf('summary_%s.csv',TAG));
assert(isfile(f), 'missing summary: %s', f);
T = readtable(f);

% The probe lives only in SSH runs: a BON run has no SSH solve to compare
% against, so its probe columns are NaN. Both curves below are therefore
% measured on the SSH trajectory, which is the point of the figure.
T = T(strcmpi(T.method,'SSH'), :);
D = T(strcmpi(T.policy,POLIC), :);

epsv = sort(unique(round(D.eps,4)), 'descend')';   % 0.25 ... 0.01
conf = 100*(1 - epsv);                             % 75 ... 99, the x axis
nE   = numel(epsv);
hBon = nan(1,nE);  hSsh = nan(1,nE);  only = nan(1,nE);  gain = nan(1,nE);
yr   = nan(1,nE);
for a = 1:nE
    r = D(abs(round(D.eps,4)-epsv(a))<1e-9, :);
    assert(height(r) > 0, 'no rows for %s eps %.2f', POLIC, epsv(a));
    yr(a)   = mean(r.hours);                       % 8760, the water year
    hBon(a) = mean(r.probe_bon_infeas);            % hours UNattained
    hSsh(a) = mean(r.probe_ssh_infeas);
    only(a) = mean(r.probe_ssh_only);
    gain(a) = mean(r.probe_obj_gain_pct);
end
assert(all(abs(yr - yr(1)) < 1e-9), ...
    'runs have different horizons (%s); the y axis is hours attained out of a fixed year', ...
    mat2str(yr));

% The probe counts hours the target is MISSED; the axis shows hours it is
% MET, so SSH sits above Bonferroni and higher reads as better. The band is
% unchanged -- complementing both curves is a reflection, not a rescaling,
% so mSsh - mBon is still exactly probe_ssh_only.
mBon = yr - hBon;
mSsh = yr - hSsh;

% ---- the three things the caption asserts, checked before anything is drawn
nCont = sum(T.probe_contain_viol);
nWors = sum(T.probe_obj_worse);
nHrs  = sum(T.hours);
assert(nCont == 0, ['CONTAINMENT BROKEN: %d hours where Bonferroni was ' ...
    'feasible and SSH was not. The caption claims zero, and without it the ' ...
    'band between the curves is no longer the recovered hours. Do not ' ...
    'publish until this is explained.'], nCont);
assert(nWors == 0, ['DOMINANCE BROKEN: %d hours where SSH returned a worse ' ...
    'objective than Bonferroni from the same state. The caption claims ' ...
    'zero.'], nWors);
assert(all(abs((mSsh - mBon) - only) < 1e-6), ...
    ['the band is not the recovered hours at every point: gap %s against ' ...
     'probe_ssh_only %s. Containment implies these are equal; if they are ' ...
     'not, the caption is wrong.'], mat2str(mSsh-mBon,5), mat2str(only,5));

% The axis is drawn as a PERCENTAGE of the water year. mBon/mSsh stay in
% hours above, so the assertions and the printed report still compare hours
% to hours; only the plotted series are rescaled.
pBon = 100*mBon./yr;
pSsh = 100*mSsh./yr;

%% ---------------------------------------------------------------- FIGURE
fig = figure('Units','inches','Position',[1 1 4 2],'Color','w');
ax  = axes(fig); hold(ax,'on');

% the band between the curves IS the hours SSH recovers -- see header
patch(ax, [conf fliplr(conf)], [pSsh fliplr(pBon)], cFill, ...
      'EdgeColor','none', 'HandleVisibility','off');

% Lines and markers are drawn as SEPARATE objects, in that order. A single
% '-o' call makes one Line object whose stroke is painted in the same pass
% as its markers, so at the exported resolution the stroke reads as sitting
% on top of the marker faces and every point looks nicked. Drawing all the
% strokes first and all the markers after puts the markers cleanly above.
% The legend is fed by NaN proxies at the end, since the real handles are
% now split in two.

plot(ax, conf, pSsh, '-', 'Color',cSSH, 'LineWidth',1.6, 'HandleVisibility','off');
plot(ax, conf, pBon, '-', 'Color',cBON, 'LineWidth',1.6, 'HandleVisibility','off');

plot(ax, conf, pSsh, 'o', 'LineStyle','none', 'MarkerSize',5, ...
     'MarkerFaceColor',cSSH, 'MarkerEdgeColor',cSSH, 'HandleVisibility','off');
plot(ax, conf, pBon, 's', 'LineStyle','none', 'MarkerSize',5, ...
     'MarkerFaceColor',cBON, 'MarkerEdgeColor',cBON, 'HandleVisibility','off');

hLeg = [plot(ax, NaN, NaN, '-o', 'Color',cSSH, 'LineWidth',1.6, ...
             'MarkerSize',5, 'MarkerFaceColor',cSSH, 'MarkerEdgeColor',cSSH),
             plot(ax, NaN, NaN, '-s', 'Color',cBON, 'LineWidth',1.6, ...
             'MarkerSize',5, 'MarkerFaceColor',cBON, 'MarkerEdgeColor',cBON)];

% NOT "feasible plan found". The lower storage bound carries a slack
% variable penalized at C.slackPenalty = 1e5 (load_case.m:27), so the LP is
% essentially never infeasible and a plan always exists. What is counted is
% hours at which the plan MEETS the joint chance constraint without invoking
% that slack -- the axis is about the reliability certificate, not about
% whether the controller produced a schedule. Keep label and caption
% consistent on this point.
ylabel(ax, {'Reliability Target Met (%)'});
xlabel(ax, 'Joint Reliability Level 1 - \epsilon (%)');
xlim(ax, [75 99]);   xticks(ax, [75 80 85 90 95 99]);
ylim(ax, [89 94]);        yticks(ax, 89:1:94);

ax.FontName = 'Times New Roman';  ax.FontSize = FS;
ax.LineWidth = 0.8;  ax.Box = 'on';  ax.TickDir = 'in';
ax.XColor = [0.30 0.30 0.30];  ax.YColor = [0.30 0.30 0.30];
ax.YGrid = 'off';  ax.XGrid = 'off';            % horizontal only, as Fig. 5
ax.GridColor = [0.85 0.85 0.85];  ax.GridAlpha = 1;
ax.Layer = 'top';                              % keep the band under the ticks

% southeast: both curves now climb to the right, so the free corner moved
legend(hLeg, {'Supporting Hyperplanes', 'Bonferroni'}, 'Location','southeast', ...
       'Box','off', 'FontSize',FS-0.5);

exportgraphics(fig, fullfile(FIGD,'ssh_vs_bon.pdf'), 'ContentType','vector');
exportgraphics(fig, fullfile(FIGD,'ssh_vs_bon.png'), 'Resolution', 400);

%% ---------------------------------------------------------------- REPORT
fprintf('%s | %s | %d SSH runs probed, %d hours\n', TAG, POLIC, height(T), nHrs);
fprintf('containment violations %d   objective-worse steps %d   (both must be 0)\n\n', ...
        nCont, nWors);
fprintf(['  1-eps    eps    BON met  SSH met   BON miss  SSH miss   recovered h' ...
         '   %% of BON miss   obj gain %%\n']);
for a = 1:nE
    fprintf('  %5.0f%%  %5.2f  %8.1f %8.1f  %9.1f %9.1f  %12.1f  %13.1f%%  %10.3f\n', ...
        conf(a), epsv(a), mBon(a), mSsh(a), hBon(a), hSsh(a), ...
        only(a), 100*only(a)/hBon(a), gain(a));
end

fprintf('\nper-scenario recovered hours (the band, scenario by scenario):\n');
for a = 1:nE
    r = sortrows(D(abs(round(D.eps,4)-epsv(a))<1e-9, :), 'scen');
    fprintf('  %3.0f%% : %s\n', conf(a), sprintf('%5d', r.probe_ssh_only));
end

k = find(abs(conf-95) < 1e-9, 1);
fprintf(['\nheadline (%s, 95%% reliability): Bonferroni meets the target in %.0f of ' ...
         'the 8760 hours,\nSSH in %.0f. The %.0f hour band is %.0f%% of the %.0f hours ' ...
         'Bonferroni misses.\nOn the hours both meet it SSH improves the objective by ' ...
         '%.2f%%, and is never worse.\n'], ...
        POLIC, mBon(k), mSsh(k), only(k), 100*only(k)/hBon(k), hBon(k), gain(k));
fprintf('wrote ssh_vs_bon.{pdf,png} to %s\n', FIGD);

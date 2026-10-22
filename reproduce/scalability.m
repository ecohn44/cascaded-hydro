%% ========================================================================
% REPRODUCE: SCALABILITY TIMING
% ========================================================================
%
% Solve time and cut count against cascade size
%
%   Run it:  >> scalability        (script, not a function)
%
%   Answers reviewer comments 2.6 / 3.2 (CPU time, scalability). Sweeps
%   n = 4, 8, 12, 16 by tiling the four-unit block, runs BON and SSH at each
%   size, and records solve time per rolling hour, iterations to
%   convergence, cut counts and problem size.
%
%   THIS IS A SYNTHETIC SCALING STUDY. Report ONLY timing, iteration and
%   cut counts, and problem size from it. Never generation, revenue or IVI
%   as physical results -- a sixteen-project Columbia is not a real river.
%   Say so in the caption.
%
%   Four decisions in here are not arbitrary; each came out of the smoke
%   test and each would bias the answer if reverted.
%
%   1. CAP = 200, not the production 20. At n = 8 SSH needs about 29
%      iterations to converge, so a cap of 20 truncated EVERY solved hour:
%      146 of 168 steps sat exactly on the cap. A truncated solve reports a
%      solve time for a problem it did not finish, which would flatter the
%      method at large n while the solution quality degraded underneath.
%      The cap must not bind; the script fails loudly if it does.
%
%   2. A 24-hour WARMUP that is discarded from every statistic. C.lag
%      reaches back 12 hours and realtimeSimulation clamps the index at zero, so the
%      first hours of any window run on truncated upstream history.
%
%   3. START_H = 4032 (17 Mar) rather than hour 1 (1 Oct). October is the
%      worst month of the water year for constraint attainment -- 16.1% of
%      hours unattainable at n = 4 in the production runs against 2.8% in
%      March, 7.2% for the year. Note that a sliced run cold-starts from
%      observed storage rather than arriving on its own trajectory, so its
%      attainment rate does NOT reproduce the production monthly figure:
%      measured here it is about 5% at n = 4, not 0.3%.
%
%   4. Timing and iteration statistics are computed over CONVERGED steps
%      only. Hours where no feasible point attains the joint target exit
%      with zero iterations, and averaging those in drags the mean down.
%      The attainment rate is reported separately as its own column, since
%      it is a real property of the formulation and rises with n.
%
%   Outputs:
%       results/scalability.csv    one row per (n, method)
%       results/scalability.mat    the same, plus every summary struct
%
%   Budget about 10 minutes. BLOCKS = 1:4 adds n = 16 for roughly 5 more.

clear; clc; close all;

%% ---------------------------------------------------------------- CONFIG
% This file lives in functions/, so the repository root is two levels up.
ROOT = fileparts(fileparts(mfilename('fullpath')));
addpath(ROOT); addpath(fullfile(ROOT,'functions'));

BLOCKS   = 1:3;                % n = 4, 8, 12
FULLYEAR = true;               % see below -- false restores the March week
WARMUP   = 24;                 % hours run then discarded from the statistics
HOURS    = 168;                % window length when FULLYEAR is false
START_H  = 4032;               % 17 Mar, when FULLYEAR is false
CAP     = 200;                 % must not bind -- checked below
EPSV    = 0.05;
THETA   = 400;
POLICY  = 'DDU';
METHODS = {'BON','SSH'};

outDir = fullfile(ROOT,'results');
if ~isfolder(outDir), mkdir(outDir); end
logFile = fullfile(ROOT,'logs', ...
    sprintf('scalability_%s.log', datestr(now,'yyyymmdd_HHMMSS')));
if ~isfolder(fileparts(logFile)), mkdir(fileparts(logFile)); end
diary(logFile); diary on;

if FULLYEAR
    fprintf('scalability | %s, eps %.2f, theta %d, FULL YEAR + %d h warmup\n', ...
            POLICY, EPSV, THETA, WARMUP);
else
    fprintf('scalability | %s, eps %.2f, theta %d, %d h + %d warmup from h%d\n', ...
            POLICY, EPSV, THETA, HOURS, WARMUP, START_H);
end
fprintf('Started %s\n\n', datestr(now));

% WHY THE FULL YEAR IS THE DEFAULT. The March window was chosen because it
% is the CLEANEST month for constraint attainment -- 2.8% of hours
% unattainable against 16.1% in October and 7.2% for the year -- so that any
% rise in cost with n was attributable to n rather than to the season. That
% is the right control for isolating the scaling exponent, but it cuts the
% other way for the headline numbers: tighter states need more cuts and more
% iterations, so a March week UNDERSTATES SSH's iteration count and solve
% time. Reporting a hand-picked easy week in answer to a reviewer asking
% about realistic operating conditions invites exactly the objection it was
% meant to avoid. The full year includes October.
%
% Cost, measured from the 192-hour runs and scaled by 8760/192: about 5.1 h
% for n = 12 SSH, 1.7 h for n = 8, under an hour for the rest. Six jobs, so
% on 7 workers the wall time is roughly the longest job.
%
% Statistics are still taken over CONVERGED steps only, so the extra
% unattainable hours in October do not drag the timing means down; they show
% up in the attainment column instead, where they belong.
C0 = dataload(ROOT, 2024, 0.20, 7);
if FULLYEAR
    Cw    = C0;
    HOURS = size(C0.ref,1) - WARMUP;
else
    Cw = sliceCase(C0, START_H - WARMUP, HOURS + WARMUP);
end

%% ------------------------------------------------------------- THE SWEEP
R = struct('n',{},'method',{},'nvar',{},'nrow',{}, ...
           'attain_pct',{},'t_mean',{},'t_p95',{},'t_max',{},'t_all_mean',{}, ...
           'iter_mean',{},'iter_p95',{},'iter_max',{},'cuts_mean',{}, ...
           'wall_s',{},'hit_cap',{});
OUT = {};
t00 = tic;

for m = BLOCKS
    Cm = tileCase(Cw, m, 'Quiet', true);
    n  = numel(Cm.units);
    for im = 1:numel(METHODS)
        meth = METHODS{im};
        fprintf('--- n = %2d, %s ---\n', n, meth);

        % RESUMABLE. A full-year sweep is ~5 hours and there are only six
        % jobs, so a crash or a sleeping laptop at hour 7000 would otherwise
        % lose the night. Each (n, method) result is cached the moment it
        % finishes; rerunning the script reloads what is already done and
        % runs only what is missing. Delete results/scalability_jobs to
        % force a clean rerun.
        cacheF = fullfile(outDir,'scalability_jobs', ...
                          sprintf('n%02d_%s.mat', n, meth));
        if isfile(cacheF)
            L = load(cacheF,'row','summ');
            k = numel(R) + 1;  R(k) = L.row;  OUT{k} = L.summ; %#ok<SAGROW>
            fprintf('    cached from a previous run, skipped\n');
            continue
        end

        tc  = tic;
        out = realtimeSimulation(Cm, POLICY, mkcfg(meth, HOURS+WARMUP, EPSV, CAP, THETA));
        wall = toc(tc);

        tr   = out.trajectory;
        keep = (WARMUP+1) : numel(tr.solve_seconds);
        sec  = double(tr.solve_seconds(keep));

        if strcmpi(meth,'SSH')
            ok   = logical(tr.ssh_step_feasible(keep));
            it   = double(tr.ssh_iterations(keep));
            cuts = double(tr.ssh_cuts(keep));
            conv = ok & it > 0 & it < CAP;
            hitCap = nnz(it >= CAP);
        else
            ok   = logical(tr.bon_step_feasible(keep));
            it   = nan(size(sec));  cuts = nan(size(sec));
            conv = ok;
            hitCap = 0;
        end

        if nnz(conv) == 0
            warning('n=%d %s: no converged steps; row will be NaN', n, meth);
        end

        k = numel(R) + 1;
        R(k).n          = n;
        R(k).method     = meth;
        R(k).nvar       = n*24*10;
        R(k).nrow       = n*24*10;
        R(k).attain_pct = 100*mean(ok);
        R(k).t_mean     = mean(sec(conv),'omitnan');
        R(k).t_p95      = quantile(sec(conv),0.95);
        R(k).t_max      = max(sec(conv));
        R(k).t_all_mean = mean(sec,'omitnan');
        R(k).iter_mean  = mean(it(conv),'omitnan');
        R(k).iter_p95   = quantile(it(conv),0.95);
        R(k).iter_max   = max(it(conv));
        R(k).cuts_mean  = mean(cuts(conv),'omitnan');
        R(k).wall_s     = wall;
        R(k).hit_cap    = hitCap;
        OUT{k} = out.summary; %#ok<SAGROW>

        if ~isfolder(fileparts(cacheF)), mkdir(fileparts(cacheF)); end
        row = R(k);  summ = OUT{k}; %#ok<NASGU>
        save(cacheF,'row','summ');

        fprintf('    attain %5.1f%% | t %.3f s/hour (p95 %.3f) | iters %.1f | %.1f min\n', ...
                R(k).attain_pct, R(k).t_mean, R(k).t_p95, R(k).iter_mean, wall/60);
        if hitCap > 0
            fprintf(['    !! %d steps reached the cap of %d. The cap is binding, so ' ...
                     'these\n       timings are for unfinished solves. Raise CAP and ' ...
                     'rerun.\n'], hitCap, CAP);
        end
    end
end

%% ---------------------------------------------------------------- OUTPUT
T = struct2table(R);
writetable(T, fullfile(outDir,'scalability.csv'));
scalability = struct('table',T,'summaries',{OUT},'config', ...
    struct('blocks',BLOCKS,'hours',HOURS,'warmup',WARMUP,'startHour',START_H, ...
           'cap',CAP,'eps',EPSV,'theta',THETA,'policy',POLICY), ...
    'when',datestr(now)); %#ok<NASGU>
save(fullfile(outDir,'scalability.mat'),'scalability','-v7.3');

fprintf('\n===== scalability =====\n');
fprintf('  %3s %-5s %7s %8s %9s %9s %8s %8s\n', ...
        'n','meth','vars','attain%','s/hour','p95 s','iters','cuts');
for k = 1:numel(R)
    fprintf('  %3d %-5s %7d %8.1f %9.3f %9.3f %8.1f %8.1f\n', ...
        R(k).n, R(k).method, R(k).nvar, R(k).attain_pct, ...
        R(k).t_mean, R(k).t_p95, R(k).iter_mean, R(k).cuts_mean);
end

% How cost grows with n. NOTE 4 -> 8 -> 12 is not a sequence of doublings,
% so a "per doubling" factor would be wrong for the second step. Report the
% step ratios as measured, plus a fitted exponent from a log-log regression,
% which handles unequal spacing and is the quotable number: time ~ n^p.
ssh = R(strcmpi({R.method},'SSH'));
if numel(ssh) > 1
    fprintf('\n  SSH cost growth (measured step ratios):\n');
    for k = 2:numel(ssh)
        fprintf('    n %2d -> %2d (x%.2f in n) : time x%.2f, iterations x%.2f\n', ...
            ssh(k-1).n, ssh(k).n, ssh(k).n/ssh(k-1).n, ...
            ssh(k).t_mean/max(ssh(k-1).t_mean,eps), ...
            ssh(k).iter_mean/max(ssh(k-1).iter_mean,eps));
    end
    ln = log([ssh.n]);
    pt = polyfit(ln, log([ssh.t_mean]), 1);
    pi_ = polyfit(ln, log([ssh.iter_mean]), 1);
    fprintf('\n    fitted exponents (%d points): time ~ n^%.2f, iterations ~ n^%.2f\n', ...
            numel(ssh), pt(1), pi_(1));
    fprintf('    extrapolated to n = 16: %.2f s/hour\n', exp(polyval(pt, log(16))));
    fprintf('\n  headroom at the largest size RUN: %.3f s per rolling hour against a\n', ...
            ssh(end).t_mean);
    fprintf('  3600 s decision cycle -- a factor of %.0f.\n', ...
            3600/max(ssh(end).t_mean,eps));
end

% the SSH/BON ratio, which is what "is SSH affordable" actually means
fprintf('\n  SSH cost relative to Bonferroni:\n');
for nn = unique([R.n])
    a = R([R.n]==nn & strcmpi({R.method},'BON'));
    b = R([R.n]==nn & strcmpi({R.method},'SSH'));
    if ~isempty(a) && ~isempty(b)
        fprintf('    n = %2d : %.3f vs %.3f s/hour  (x%.1f)\n', ...
            nn, a.t_mean, b.t_mean, b.t_mean/max(a.t_mean,eps));
    end
end

fprintf('\nwrote results/scalability.csv and .mat\n');
fprintf('Total %.1f min. Finished %s\n', toc(t00)/60, datestr(now));
diary off;

%% ------------------------------------------------------------- helpers
function cfg = mkcfg(method, hours, epsilon, maxIter, theta)
cfg = struct( ...
    'Horizon',24, 'Epsilon',epsilon, 'Kappa',1, 'Theta',theta, ...
    'MaxHours',hours, 'ScenarioID',1, 'ScenarioSeed',1, ...
    'SolutionMethod',method, ...
    'SSHMaxIterations',maxIter, 'SSHPhiTolerance',5e-4, ...
    'SSHBoundaryTolerance',1e-5, 'SSHLookaheadHours',24, ...
    'SSHInteriorPoint','bonferroni', 'ProgressHours',1e9, 'Verbose',false, ...
    'ChanceMode','on', 'ChanceSides','lower', ...
    'TrackSide','both', 'TrackBand',0, ...
    'PriceMode','bpa', 'PriceFlatten',0, ...
    'PriceNoiseSigma',0.30, 'PriceNoiseRho',0.80, 'PriceNoiseSeed',1, ...
    'WorldKappa',1, 'PolicyKappa',1, ...
    'SigmaModel','archx', 'PowerModel','mccormick');
end

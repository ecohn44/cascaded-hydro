%% ========================================================================
% TEST: FULL-YEAR SMOKE
% ========================================================================
%
% Two checks before committing five hours overnight
%
%   Run it:  >> smokeTestFullYear        (script, about 8 minutes)
%
%   reproduce/scalability with FULLYEAR = true is NOT resumable -- there are no
%   per-job files, so a failure at hour 7000 loses the whole run. These are
%   the two things that differ from the 168-hour March run and have never
%   been exercised.
%
%   CHECK 1: tileCase on an unsliced year. Every previous call tiled a
%   sliced window. tileCase does not reference sliceCase's fields, so this
%   should be fine, but "should be fine" is what an overnight run is for
%   finding out about. Instant.
%
%   CHECK 2: the iteration cap in the WORST month. This is the real risk.
%   The cap exists because a truncated solve reports a time for a problem it
%   did not finish, which flatters the method while quality degrades
%   underneath. In March at n = 12 the maximum was 68 iterations against a
%   cap of 200 -- comfortable. But March is the cleanest month for
%   attainment (2.8% unattainable) and October the worst (16.1%), and
%   tighter states need more cuts. The full year includes October. If the
%   cap binds there, every timing number in the run is for unfinished
%   solves and the whole night is wasted.
%
%   The window is h504 (22 Oct), which the monthly scan in sliceCase
%   identifies as the WORST 168 hours of the water year at 24.7%
%   unattainable. If the cap survives that, it survives the year.
%
%   Also reports measured wall time so the overnight estimate can be
%   checked against the hard month rather than the easy one.

clear; clc;

% This test lives in tests/, so the repository root is one level up.
ROOT = fileparts(fileparts(mfilename('fullpath')));
addpath(ROOT);
addpath(fullfile(ROOT,'functions'));

BLOCK   = 3;        % n = 12, the largest and slowest size
WORST_H = 504;      % 22 Oct, the worst 168-hour window of WY2024
HOURS   = 168;
WARMUP  = 24;
CAP     = 200;      % same as reproduce/scalability

C0 = dataload(ROOT, 2024, 0.20, 7);
T  = size(C0.ref,1);

%% ------------------------------------------------- CHECK 1: tile the year
fprintf('check 1: tileCase on the unsliced year\n');
tc = tic;
Cfull = tileCase(C0, BLOCK, 'Quiet', true);
assert(numel(Cfull.units) == 4*BLOCK, ...
    'tiled case has %d units, expected %d', numel(Cfull.units), 4*BLOCK);
assert(size(Cfull.ref,1) == T, ...
    'tiled case has %d hours, expected %d -- tiling changed the row count', ...
    size(Cfull.ref,1), T);
assert(all(isfinite(Cfull.qLow(:))) && all(isfinite(Cfull.qHigh(:))), ...
    'tiled envelope contains non-finite values');
fprintf('  PASS: %d units, %d hours, envelope finite (%.1f s)\n\n', ...
        numel(Cfull.units), size(Cfull.ref,1), toc(tc));
clear Cfull

%% --------------------------------------- CHECK 2: the cap in the worst week
fprintf('check 2: n = %d SSH on the worst week of the year (h%d, 22 Oct)\n', ...
        4*BLOCK, WORST_H);
Cw = sliceCase(C0, WORST_H - WARMUP, HOURS + WARMUP);
Cm = tileCase(Cw, BLOCK, 'Quiet', true);

tc  = tic;
out = realtimeSimulation(Cm, 'DDU', mkcfg('SSH', HOURS+WARMUP, 0.05, CAP, 400));
wall = toc(tc);

tr   = out.trajectory;
keep = (WARMUP+1) : numel(tr.solve_seconds);
it   = double(tr.ssh_iterations(keep));
ok   = logical(tr.ssh_step_feasible(keep));
sec  = double(tr.solve_seconds(keep));
conv = ok & it > 0 & it < CAP;
hitCap = nnz(it >= CAP);

fprintf('  attainment   %.1f%%  (March was 82.7%%)\n', 100*mean(ok));
fprintf('  iterations   mean %.1f, p95 %.1f, max %d   (March max was 68)\n', ...
        mean(it(conv)), quantile(it(conv),0.95), max(it(conv)));
fprintf('  solve time   mean %.3f s, p95 %.3f s       (March mean was 2.400)\n', ...
        mean(sec(conv)), quantile(sec(conv),0.95));
fprintf('  wall         %.1f min for %d hours\n', wall/60, HOURS+WARMUP);
fprintf('  extrapolated full year at n = %d: %.1f h\n\n', ...
        4*BLOCK, wall*(size(C0.ref,1)/(HOURS+WARMUP))/3600);

if hitCap > 0
    error('smokeTestFullYear:capBinds', ...
        ['FAIL: %d of %d steps reached the cap of %d in the worst week. ' ...
         'Every timing number in an overnight run would be for unfinished ' ...
         'solves. Raise CAP in reproduce/scalability.m to at least %d and rerun ' ...
         'this test before starting the full year.'], ...
        hitCap, numel(it), CAP, 2*max(it));
end

fprintf(['PASS: the cap does not bind in the worst week (max %d of %d).\n' ...
         'reproduce/scalability with FULLYEAR = true is safe to start.\n'], ...
        max(it(conv)), CAP);

%% ------------------------------------------------------------- helper
function cfg = mkcfg(method, hours, epsilon, maxIter, theta)
% copied from reproduce/scalability.m -- keep the two in step or this tests
% something other than what runs overnight
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

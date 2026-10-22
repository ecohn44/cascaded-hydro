%% ========================================================================
% TEST: TILED CASCADE SMOKE
% ========================================================================
%
% Check the tiled cascade before timing anything with it
%
%   Run it:  >> smokeTestScale        (script, not a function)
%
%   A. IDENTITY. tileCase(C,1) through realtimeSimulation must equal C
%      through realtimeSimulation, bit for bit.
%   B/C/D. n = 8 at the headline epsilon, at a raised iteration cap, and at
%      a matched per-unit tolerance.
%
%   READ THIS BEFORE TRUSTING ANY COLUMN NAMED "unattainable".
%   realtimeSimulation.m line 854 defines
%       SSHTargetUnattainable(t) = strcmp(solutionMethod,'SSH') && ~SSHStepFeasible(t)
%   and for an SSH run chance_infeasible_steps is nnz(~SSHStepFeasible). The
%   two are therefore THE SAME NUMBER by construction, not two independent
%   measurements that happen to agree. The field name promises a
%   distinction -- "no feasible point reaches the target" -- that it does
%   not actually make. An earlier version of this script compared them and
%   concluded "unattainable"; that conclusion was vacuous.
%
%   What actually distinguishes the failure modes is WHERE in
%   applySSHNative a step gives up, and the iteration count is the tell:
%       0 iterations, feasible     Bonferroni already met the target, no
%                                  cuts required (early return, line ~52)
%       0 iterations, infeasible   the base SSH LP was not OPTIMAL
%                                  (line ~43)
%       at the cap, infeasible     ran the full budget and never converged
%                                  (line ~208, falls back)
%       below the cap, feasible    converged normally
%   So the diagnostic below cross-tabulates iterations against
%   ssh_step_feasible instead of comparing two names for one quantity.
%
%   It also saves every run to results/tests/scale_smoke.mat, so this can
%   be re-asked without re-running anything.
%
%   Counts must NOT be compared across B, C and D. Each is a separate
%   closed-loop simulation: change the cap or epsilon and the trajectory
%   diverges, so a different set of hours is visited. Compare rates and
%   distributions within a run, not hour counts between runs.

clear; clc; close all;

%% ---------------------------------------------------------------- CONFIG
% This test lives in tests/, so the repository root is one level up.
ROOT = fileparts(fileparts(mfilename('fullpath')));
addpath(ROOT);
addpath(fullfile(ROOT,'functions'));

YEAR = 2024;  FLOORPCT = 0.20;  HALFWIN = 7;
HOURS   = 168;
METHOD  = 'SSH';
POLICY  = 'DDU';
ITER_HI = 60;
START_H = 4032;    % 17 Mar. Hour 1 is 1 Oct, the worst month of the year for
                   % constraint attainment (16.1% of hours unattainable at n=4
                   % in the production runs, against 2.8% in March and 7.2% for
                   % the year). See sliceCase.m for the full monthly table.
                   % Set START_H = 1 to reproduce the October behaviour.

fprintf('smokeTestScale  |  %s %s, %d h from hour %d, theta 400\n\n', ...
        POLICY, METHOD, HOURS, START_H);
C  = dataload(ROOT, YEAR, FLOORPCT, HALFWIN);
C  = sliceCase(C, START_H, HOURS);     % slice before tiling; they commute
nb = numel(C.units);

%% ================================================== A. IDENTITY AT m = 1
fprintf('--- A. identity ---\n');
C1 = tileCase(C, 1, 'Quiet', true);
ref = realtimeSimulation(C,        POLICY, mkcfg(METHOD,HOURS,0.05,20));
new = realtimeSimulation(C1, POLICY, mkcfg(METHOD,HOURS,0.05,20));
flds = {'energy_MWh','lower_ivi_kaf','spill_kaf','curtail_MWh', ...
        'lower_exit_pct','joint_reliability_mean'};
passA = true;
for k = 1:numel(flds)
    a = ref.summary.(flds{k}); b = new.summary.(flds{k});
    if abs(a-b)/max(1,abs(a)) > 1e-12, passA = false; end
end
dV = max(abs(ref.trajectory.V(:) - new.trajectory.V(:)));
passA = passA && dV <= 1e-9;
fprintf('  six metrics + storage trajectory: %s (max |dV| %.1e)\n\n', tf(passA), dV);
if ~passA, fprintf('STOP. The n=4 path changed.\n'); return, end

%% ============================================ B/C/D. the n = 8 questions
C8 = tileCase(C, 2, 'Quiet', true);
fprintf('--- B. n=8, eps 0.05, cap 20 ---\n');
r8 = realtimeSimulation(C8, POLICY, mkcfg(METHOD,HOURS,0.05,20));
fprintf('--- C. n=8, eps 0.05, cap %d ---\n', ITER_HI);
rC = realtimeSimulation(C8, POLICY, mkcfg(METHOD,HOURS,0.05,ITER_HI));
fprintf('--- D. n=8, eps 0.10, cap 20 ---\n');
rD = realtimeSimulation(C8, POLICY, mkcfg(METHOD,HOURS,0.10,20));

cases = {'n=4  e.05 c20', ref, 20; ...
         'n=8  e.05 c20', r8,  20; ...
         'n=8  e.05 c60', rC,  ITER_HI; ...
         'n=8  e.10 c20', rD,  20};

fprintf('\n  %-14s %9s %8s %8s %8s %8s\n', ...
        'case','IVI kaf','exit %','s/hour','p95 s','iters');
for k = 1:size(cases,1)
    s = cases{k,2}.summary;
    fprintf('  %-14s %9.2f %8.3f %8.3f %8.3f %8.2f\n', cases{k,1}, ...
        s.lower_ivi_kaf, s.lower_exit_pct, s.solve_time_mean_s, ...
        s.solve_time_p95_s, s.ssh_iterations_mean);
end

%% ------------------------------------- iterations x feasibility crosstab
fprintf('\n  where the SSH steps land (hours out of %d):\n', HOURS);
fprintf('  %-14s %10s %10s %12s %12s %10s\n', 'case', ...
        '0 iter ok', '0 iter BAD', '<cap conv', 'AT CAP bad', 'AT CAP ok');
for k = 1:size(cases,1)
    tr  = cases{k,2}.trajectory;
    cap = cases{k,3};
    it  = double(tr.ssh_iterations(:));
    ok  = logical(tr.ssh_step_feasible(:));
    fprintf('  %-14s %10d %10d %12d %12d %10d\n', cases{k,1}, ...
        nnz(it==0 &  ok), nnz(it==0 & ~ok), ...
        nnz(it>0 & it<cap &  ok), nnz(it>=cap & ~ok), nnz(it>=cap &  ok));
end

% iterations needed by the steps that actually converged -- the quantity
% that scales, and the one a capped run cannot report
fprintf('\n  iterations among converged steps (mean / p95 / max):\n');
for k = 1:size(cases,1)
    tr  = cases{k,2}.trajectory;
    cap = cases{k,3};
    it  = double(tr.ssh_iterations(:));
    ok  = logical(tr.ssh_step_feasible(:));
    sel = it>0 & ok & it<cap;
    if nnz(sel)==0
        fprintf('  %-14s  none converged below the cap\n', cases{k,1});
    else
        v = it(sel);
        fprintf('  %-14s  %6.1f / %4d / %4d   (n = %d)\n', cases{k,1}, ...
            mean(v), round(quantile(v,0.95)), max(v), numel(v));
    end
end

%% ---------------------------------------------------------------- SAVE
outDir = fullfile(ROOT,'results','tests');   % test output, not a result
if ~isfolder(outDir), mkdir(outDir); end
scale_smoke = struct('ref',ref,'r8',r8,'rC',rC,'rD',rD, ...
    'hours',HOURS,'iterHi',ITER_HI,'when',datestr(now)); %#ok<NASGU>
save(fullfile(outDir,'scale_smoke.mat'),'scale_smoke','-v7.3');
fprintf('\nsaved results/tests/scale_smoke.mat\n');

%% ------------------------------------------------------------ SANITY
thr4 = sum(ref.trajectory.U,'all')*C.KAF;
thr8 = sum(r8.trajectory.U, 'all')*C.KAF;
fprintf('spill as a share of throughput: %.4f%% at n=4, %.4f%% at n=8\n', ...
        100*ref.summary.spill_kaf/max(thr4,eps), ...
        100*r8.summary.spill_kaf /max(thr8,eps));
for b = 1:C8.tiledBlocks
    cols = (b-1)*nb + (1:nb);
    fprintf('  block %d  spill %7.2f kaf   mean V %s\n', b, ...
        sum(r8.trajectory.spill(:,cols),'all')*C.KAF, ...
        sprintf('%7.1f', mean(r8.trajectory.V(:,cols),1)));
end

%% ------------------------------------------------------------- helpers
function cfg = mkcfg(method, hours, epsilon, maxIter)
cfg = struct( ...
    'Horizon',24, 'Epsilon',epsilon, 'Kappa',1, 'Theta',400, ...
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

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end

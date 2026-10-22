%% ========================================================================
% REPRODUCE: MPC BENCHMARKS
% ========================================================================
%
% Deterministic and robust MPC at the headline operating point
%
%   Run it:  >> benchmarks        (script, not a function)
%
%   Answers reviewer comment 1.8: "the proposed approach should be
%   quantitatively compared against well-established real-time dispatch
%   strategies ... including deterministic MPC, scenario-based stochastic
%   MPC, and robust MPC formulations." All three are here.
%
%   SCENARIO MPC IS NOT EXPENSIVE, contrary to the first draft of
%   claude/PLAN_mpc_benchmarks.md, which claimed the LP grows linearly in S
%   and budgeted 3-7 h per full-year run. It does not grow at all. Each
%   sample contributes a separate storage row, so enforcing all S samples at
%   all (i,k) is the same as enforcing the worst sample at each (i,k); and
%   since the margin is z*sqrt(k)*sigma while the shock is
%   sqrt(k)*sigma*zeta, sigma cancels and the whole scenario program reduces
%   to choosing z as an order statistic of the standardised shock. See
%   scenarioMargin.m for the three conditions that make the collapse valid.
%   Scenario MPC therefore costs the same per job as every other row here.
%
%   RUN verifyParity FIRST. It is the gate that realtimeSimulation
%   still reproduces realtimeSimulation. Without it these rows may be comparing our
%   controller against a subtly different version of itself.
%
%   ---------------------------------------------------------------------
%   SHOULD THE BENCHMARKS USE DIU OR DDU? Both, and the answer matters.
%
%   The two axes are orthogonal. DIU-vs-DDU is what sigma the controller
%   BELIEVES; the margin rule is how that sigma ENTERS the constraint. The
%   reviewer asked about the second axis.
%
%   Run DIU alone and the benchmark is faithful to practice -- no operator
%   models sigma as a function of their own release ramp, so DIU is what
%   "widely adopted in literature and industry" actually means. But then our
%   method differs from every baseline on BOTH axes at once, and a reviewer
%   can say the whole gain is the better sigma rather than the method. That
%   objection would be fair and we would have no answer.
%
%   So every margin rule that uses sigma is run in both, and the table reads
%   as a grid:
%       DIU -> DDU at fixed margin rule  = the decision-dependence gain
%       BON -> SSH at fixed sigma model  = the allocation gain
%   The DIU rows are the honest baselines to put in the paper; the DDU rows
%   are what separate the two contributions. Quote the DIU rows in the text
%   and keep the DDU rows for the reviewer who asks.
%
%   DETERMINISTIC MPC IS THE EXCEPTION and gets one row, not two. At z = 0
%   sigma never enters the constraint, so DIU and DDU are the same
%   controller -- already verified on the old benchmark, identical to 0 MWh
%   in all 7 scenarios. Running both would be two names for one run.
%
%   THE CHANCE-CONSTRAINED ROWS ARE NOT RE-RUN. BON and SSH, DIU and DDU, at
%   theta = 400 and eps = 0.05 already exist in BPA_HL2_TH0000400. Pull them
%   from that summary when assembling the table. Re-running them here would
%   only risk a different configuration.
%   ---------------------------------------------------------------------
%
%   WHY GAMMA IS NOT DERIVED FROM EPSILON. The robust margin is
%   Gamma*sqrt(k)*sigma, algebraically the same row as Bonferroni with z
%   replaced by Gamma -- which is exactly the point. Gamma is the size of
%   the uncertainty set and carries no probabilistic calibration, so the
%   realised reliability it delivers is whatever the hydrology hands back.
%   The grid brackets the two multipliers the paper already uses:
%       eps = 0.05, 4 lower tails -> z = 2.2414
%       eps = 0.01, 4 lower tails -> z = 2.8070
%   so Gamma in [1, 3] spans from clearly under-protective to beyond our
%   most conservative chance-constrained cell. Do NOT set Gamma from eps;
%   that collapses robust MPC onto Bonferroni and destroys the comparison.
%
%   EVERY PHYSICS SETTING IS INHERITED from the headline sweep's own
%   sweep_config.mat rather than retyped here, and the shocks come from the
%   same scenario bank and indices, so these rows are paired with the
%   headline rows scenario by scenario. Changing a setting in this file
%   instead of inheriting it is how a benchmark stops being like-for-like.
%
%   THE SCENARIO GRID, and why both variants are run. With NDiscard = 0 the
%   multiplier is the sample MAXIMUM, which is the pure scenario program: it
%   grows without bound in S, because it enforces every drawn sample and the
%   Gaussian support is unbounded. More samples make it WORSE, which is the
%   opposite of the usual intuition and is the finding worth reporting. With
%   NDiscard = floor(eps*S) it is sampling-with-discarding (Campi-Garatti),
%   which targets the true joint quantile and is the honest competitor to
%   SSH. Reporting only the discarding variant would hide the first result;
%   reporting only the pure variant would be a straw man. Run both.
%
%   Before committing any compute, run `scenarioMargin` with no arguments.
%   It prints the multiplier each (S, NDiscard) chooses in a few seconds,
%   with no dispatch solves, next to the Bonferroni z for reference. If the
%   discarding column has not settled by S = 1000, extend the grid there
%   rather than discovering it after a 70 minute sweep.
%
%   Resumable: an existing job file is skipped, so a crashed run resumes.
%
%   Budget: 92 jobs at roughly 5 min each -- about 7.7 h serial, closer to
%   70 min on the 7 workers the headline sweep used.
%
%   Outputs:
%       results/runs/BENCH_TH0000400_E050/jobs/*.mat
%       results/runs/benchmarks.csv

clear; clc;

%% ---------------------------------------------------------------- CONFIG
% This file lives in functions/, so the repository root is two levels up.
ROOT = fileparts(fileparts(mfilename('fullpath')));
addpath(ROOT); addpath(fullfile(ROOT,'functions'));

SRC_TAG = 'BPA_HL2_TH0000400';          % inherit physics from here
TAG     = 'BENCH_TH0000400_E050';
EPSV    = 0.05;
SCEN    = 1:4;
% Robust uncertainty-set sizes. The sweep spans under- to over-protective,
% and 2.2340 is included because it is the Gamma whose implied per-step
% joint reliability is exactly the requested 1 - eps:
%     Phi(Gamma)^n = 0.95  ->  Gamma = Phi^-1(0.95^(1/4)) = 2.2340
% That is the row the benchmark table reports, so robust MPC is compared at
% the same reliability as the other formulations rather than at whatever a
% round Gamma happens to deliver.
%
% SAY IN THE TEXT WHAT THIS COSTS. Reaching 2.2340 requires computing the
% probability the set implies, which is the chance-constrained calculation.
% Robust optimization has no route to it on its own, so the calibrated row
% is evidence that robust MPC on this problem IS the chance constraint,
% not evidence that it can be calibrated in practice. Report the implied
% reliability of a hand-picked Gamma alongside it (Gamma = 2 gives 0.912)
% or the table quietly grants robust MPC a capability it does not have.
GAMMAS  = [1.0 1.5 2.0 2.2340 2.5 3.0];
NSCENS  = [20 100 1000];              % scenario-program sample counts
WORKERS = 7;

outDir = fullfile(ROOT,'results','runs',TAG);
jobDir = fullfile(outDir,'jobs');
if ~isfolder(jobDir), mkdir(jobDir); end

%% ------------------------------------------------- inherit the headline cfg
cfgFile = fullfile(ROOT,'results','runs',SRC_TAG,'sweep_config.mat');
assert(isfile(cfgFile), ['missing %s -- the benchmarks inherit every physics ' ...
    'setting from the headline sweep so the comparison stays like-for-like'], cfgFile);
S = load(cfgFile,'cfg');  cfg = S.cfg;
assert(abs(cfg.Theta-400) < 1e-9, ...
    'inherited theta is %g, expected 400; wrong source tag?', cfg.Theta);

fprintf('%s | inheriting physics from %s\n', TAG, SRC_TAG);
fprintf('  theta %g, year %d, floor P%02.0f, price %s, sigma %s, power %s\n', ...
        cfg.Theta, cfg.Year, 100*cfg.FloorPct, cfg.PriceMode, cfg.SigmaModel, ...
        cfg.PowerModel);
fprintf('  eps %.2f, scenarios %s\n', EPSV, mat2str(SCEN));
fprintf('  z for reference: eps 0.05 -> %.4f, eps 0.01 -> %.4f\n\n', ...
        -sqrt(2)*erfcinv(2*(1-0.05/4)), -sqrt(2)*erfcinv(2*(1-0.01/4)));

% The stored config records an ABSOLUTE BankFile path, captured on the machine
% that produced the tag. It therefore resolves nowhere else, and no longer
% resolves here either now that the shock bank lives in data/. Resolve it from
% ROOT instead; the bank itself is fixed, so this changes no number.
cfg.BankFile = fullfile(ROOT,'data','scenario_bank_WY2024_N10_seed20260920_T8760.mat');
assert(isfile(cfg.BankFile), 'missing the shock bank: %s', cfg.BankFile);

C0 = dataload(ROOT, cfg.Year, cfg.FloorPct, 7);
B  = load(cfg.BankFile,'bank');  bank = B.bank;

%% ------------------------------------------------------------- the job list
jobs = struct('scen',{},'rule',{},'val',{},'ndisc',{},'policy',{},'kappa',{},'file',{});
for s = SCEN
    % deterministic MPC: z = 0, so sigma never enters and DIU == DDU.
    % One row, not two.
    jobs(end+1) = mkjob(s,'none',0,0,'DIU',0,jobDir); %#ok<SAGROW>

    % robust: both sigma models, for the reason in the header
    for g = GAMMAS
        jobs(end+1) = mkjob(s,'robust',g,0,'DIU',0,jobDir); %#ok<SAGROW>
        jobs(end+1) = mkjob(s,'robust',g,0,'DDU',1,jobDir); %#ok<SAGROW>
    end

    % scenario program: pure (r = 0) and with discarding (r = floor(eps*S)).
    % unique() guards the small-S case where floor(eps*S) is itself 0 and
    % the two variants coincide -- running it twice would just be a slower
    % way of getting the same number.
    for S = NSCENS
        for r = unique([0 floor(EPSV*S)])
            jobs(end+1) = mkjob(s,'scenario',S,r,'DIU',0,jobDir); %#ok<SAGROW>
            jobs(end+1) = mkjob(s,'scenario',S,r,'DDU',1,jobDir); %#ok<SAGROW>
        end
    end
end

done    = arrayfun(@(j) isfile(j.file), jobs);
pending = jobs(~done);
fprintf('%d jobs, %d done, %d pending\n\n', numel(jobs), nnz(done), numel(pending));

%% ------------------------------------------------------------------- run
t00 = tic;
if ~isempty(pending)
    if WORKERS > 1
        if isempty(gcp('nocreate')), parpool(WORKERS); end
        parfor q = 1:numel(pending)
            bench_run_one(pending(q), C0, bank, cfg, EPSV);
        end
    else
        for q = 1:numel(pending)
            bench_run_one(pending(q), C0, bank, cfg, EPSV);
        end
    end
end
fprintf('\nall jobs finished in %.1f min\n', toc(t00)/60);

%% ---------------------------------------------------------------- collect
rows = {};
for k = 1:numel(jobs)
    if ~isfile(jobs(k).file), continue; end
    L = load(jobs(k).file,'summary','job');
    sm = L.summary;  jb = L.job;
    rows{end+1} = { jb.scen, jb.rule, jb.val, jb.ndisc, jb.policy, jb.kappa, ...
        getf(sm,'margin_z'), getf(sm,'energy_MWh'), getf(sm,'lower_ivi_kaf'), ...
        getf(sm,'lower_exit_pct'), getf(sm,'chance_infeasible_steps'), ...
        getf(sm,'planned_slack_kaf'), getf(sm,'spill_kaf'), ...
        getf(sm,'solve_time_mean_s'), getf(sm,'solve_time_p95_s') }; %#ok<SAGROW>
end
assert(~isempty(rows), 'no finished jobs to collect');
T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'scen','rule','val','ndisc','policy','kappa','margin_z','energy_MWh', ...
     'lower_ivi_kaf','lower_exit_pct','unattained_steps','planned_slack_kaf', ...
     'spill_kaf','solve_mean_s','solve_p95_s'});
writetable(T, fullfile(ROOT,'results','runs','benchmarks.csv'));

%% ---------------------------------------------------------------- report
fprintf('\n===== benchmarks, theta 400, eps %.2f, mean over %d scenarios =====\n', ...
        EPSV, numel(SCEN));
fprintf('  %-9s %8s %6s %-5s %9s %12s %11s %9s %11s\n', ...
        'rule','G or S','discard','sigma','margin z','energy GWh','IVI kaf', ...
        'exit %','unattained');
key = strcat(T.rule, '|', string(T.val), '|', string(T.ndisc), '|', T.policy);
for u = unique(key,'stable')'
    g = T(key==u, :);
    fprintf('  %-9s %8g %6d %-5s %9.4f %12.1f %11.0f %9.3f %11.1f\n', ...
        g.rule{1}, g.val(1), g.ndisc(1), g.policy{1}, mean(g.margin_z), ...
        mean(g.energy_MWh)/1e3, mean(g.lower_ivi_kaf), ...
        mean(g.lower_exit_pct), mean(g.unattained_steps));
end
fprintf(['\nThe margin z column is the whole benchmark in one number: it is what\n' ...
         'each method chooses in place of the Bonferroni %.4f, and the exit %%\n' ...
         'column is what that choice actually delivers.\n'], ...
        -sqrt(2)*erfcinv(2*(1-EPSV/4)));
fprintf(['\nThe chance-constrained rows (BON/SSH x DIU/DDU) are NOT here: they\n' ...
         'already exist in %s at this exact cell. Pull them from that summary\n' ...
         'when assembling the paper table.\n'], SRC_TAG);
fprintf('wrote results/runs/benchmarks.csv\n');

%% ------------------------------------------------------------- helpers
function j = mkjob(scen, rule, val, ndisc, policy, kappa, jobDir)
% val is Gamma for 'robust', the sample count S for 'scenario', unused for
% 'none'. The filename has to distinguish every arm or a resumed run will
% silently skip a job it never did.
switch rule
    case 'robust',   tok = sprintf('g%03d', round(100*val));
    case 'scenario', tok = sprintf('S%05d_r%04d', val, ndisc);
    otherwise,       tok = 'det';
end
f = sprintf('s%02d_%s_%s_%s.mat', scen, rule, tok, policy);
j = struct('scen',scen,'rule',rule,'val',val,'ndisc',ndisc, ...
           'policy',policy,'kappa',kappa,'file',fullfile(jobDir,f));
end

function bench_run_one(job, C0, bank, cfg, epsv)
t0 = tic;
try
    C = C0;  C.zseq = bank.zseq(:,:,job.scen);

    % Copied from sweepRunOne.m so the benchmark rows are built exactly
    % like the headline rows. Two deliberate departures, both asserted:
    % SolutionMethod is forced to BON, and MarginRule/Gamma are added.
    cr = struct('Horizon',cfg.Horizon,'Epsilon',epsv,'Kappa',job.kappa, ...
        'Theta',cfg.Theta,'MaxHours',cfg.MaxHours,'ScenarioID',job.scen, ...
        'ScenarioSeed',bank.seeds(job.scen),'SolutionMethod','BON', ...
        'SSHMaxIterations',20,'SSHPhiTolerance',5e-4, ...
        'SSHBoundaryTolerance',1e-5,'SSHLookaheadHours',cfg.Horizon, ...
        'SSHInteriorPoint','bonferroni','ProgressHours',1e9, ...
        'Verbose',false,'ChanceSides',cfg.ChanceSides,'TrackBand',cfg.TrackBand, ...
        'TrackSide',cfg.TrackSide,'PriceMode',cfg.PriceMode, ...
        'SigmaModel',cfg.SigmaModel,'PriceFlatten',cfg.PriceFlatten, ...
        'PriceNoiseSigma',cfg.PriceNoiseSigma,'PriceNoiseRho',cfg.PriceNoiseRho, ...
        'PriceNoiseSeed',cfg.PriceNoiseSeed,'ChanceMode',cfg.ChanceMode, ...
        'WorldKappa',cfg.WorldKappa,'PolicyKappa',cfg.PolicyKappa, ...
        'PowerModel',cfg.PowerModel, ...
        'MarginRule',job.rule);
    switch job.rule
        case 'robust'
            cr.Gamma = job.val;
        case 'scenario'
            cr.NScen    = job.val;
            cr.NDiscard = job.ndisc;
            % fixed across inflow scenarios on purpose: the multiplier is a
            % property of the sampling design, not of the hydrology, so
            % every scenario in an arm must draw the identical z
            cr.ScenSeed = 20260928;
    end

    % SSH would drive the solution back to the 1-eps target with cuts,
    % overriding whatever margin rule we asked for. BON is the only
    % solution method for which the margin rule IS the constraint.
    assert(strcmp(cr.SolutionMethod,'BON'), ...
        'benchmark rows must use BON; SSH cuts would override the margin rule');

    out = realtimeSimulation(C, job.policy, cr);
    summary = out.summary; %#ok<NASGU>
    j = job; %#ok<NASGU>
    save_bench(job.file, summary, j);
    fprintf('[%s] DONE s%d %s %g/%d %s: %.0f s\n', datestr(now,'HH:MM'), ...
            job.scen, job.rule, job.val, job.ndisc, job.policy, toc(t0));
catch ME
    fid = fopen(strrep(job.file,'.mat','_error.txt'),'w');
    fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off'));  fclose(fid);
    fprintf('[%s] FAILED s%d %s %g/%d %s: %s\n', datestr(now,'HH:MM'), ...
            job.scen, job.rule, job.val, job.ndisc, job.policy, ME.message);
end
end

function save_bench(f, summary, job) %#ok<INUSD>
save(f,'summary','job');
end

function v = getf(s, f)
if isfield(s,f), v = double(s.(f)); else, v = NaN; end
end

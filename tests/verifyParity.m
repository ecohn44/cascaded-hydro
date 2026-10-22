%% ========================================================================
% TEST: PARITY GATE
% ========================================================================
%
% Does realtimeSimulation still reproduce the published run?
%
%   Run it:  >> verifyParity        (script, not a function)
%
%   THIS IS A GATE, NOT A TEST. realtimeSimulation.m carries the MarginRule
%   knob that the benchmark table varies. At cfg.MarginRule = 'chance' -- the
%   default, and the paper's controller -- it must reproduce the published run
%   BIT FOR BIT. If it does not, the knob is not inert, and every benchmark
%   number produced from it compares the controller against a subtly different
%   version of itself rather than against deterministic or robust MPC. Nothing
%   downstream is trustworthy until this prints PASS.
%
%   HOW IT AVOIDS FOOLING ITSELF. The configuration is not retyped here. It
%   is read back from the headline sweep's own sweep_config.mat and rebuilt
%   with the identical expression sweepRunOne.m uses, so a default that
%   changed in runSweep.m since the headline ran would show up as a
%   mismatch instead of being silently reproduced. The shock path comes from
%   the same scenario bank and the same scenario index, so the realised world
%   is identical too.
%
%   WHICH JOB. s01, eps = 0.05, DDU kappa = 1, BON -- about five minutes.
%   BON rather than SSH deliberately: SSH is a 27 minute run and adds an
%   iterative solver whose tolerances make exact equality a weaker test than
%   it looks. BON exercises the whole margin path, which is the only thing
%   the patch touches.
%
%   WHAT IS COMPARED. Every numeric scalar field the two summaries share,
%   not a hand-picked few, so a field nobody thought about still fails the
%   gate. Timing fields are excluded because they are not reproducible.
%
%   Output: prints a per-field report and errors on the first mismatch.

clear; clc;

% This test lives in tests/, so the repository root is one level up.
ROOT = fileparts(fileparts(mfilename('fullpath')));
addpath(ROOT);
addpath(fullfile(ROOT,'functions'));

TAG  = 'BPA_HL2_TH0000400';
JOB  = 's01_e050_DDU_k100_BON.mat';

% Timing is wall-clock and never reproduces; everything else must.
%
% MATCHED BY PATTERN, NOT BY NAME. The first version of this gate listed six
% timing fields explicitly and missed solve_time_median_s, which then failed
% the gate at 5e-3 and looked like a real defect for as long as it took to
% read the field name. An enumerated list silently turns every field someone
% adds later into a spurious failure, which is the failure mode that trains
% people to ignore the gate. Anything whose name mentions solve time, wall
% time or seconds is excluded; everything else must match exactly.
SKIP_PAT = 'solve|wall|elapsed|seconds|_time_|time_s$';
TOL  = 1e-9;          % relative, on fields whose magnitude is above 1

%% ------------------------------------------------------------- the stored run
tagDir = fullfile(ROOT,'results','runs',TAG);
jobFile = fullfile(ROOT,'tests','reference',JOB);   % shipped reference run
cfgFile = fullfile(tagDir,'sweep_config.mat');
assert(isfile(jobFile), 'missing job: %s', jobFile);
assert(isfile(cfgFile), 'missing sweep config: %s', cfgFile);

J = load(jobFile,'summary','job');
S = load(cfgFile,'cfg');
cfg = S.cfg;  job = J.job;  ref = J.summary;

fprintf('parity gate\n');
fprintf('  reference : %s / jobs/%s\n', TAG, JOB);
fprintf('  job       : scen %d, eps %.2f, %s kappa %.2f, %s\n\n', ...
        job.scen, job.eps, job.policy, job.kappa, job.method);

%% ------------------------------------------------------------- rebuild the run
% The stored config records an ABSOLUTE BankFile path, captured on the machine
% that produced the tag. It therefore resolves nowhere else, and no longer
% resolves here either now that the shock bank lives in data/. Resolve it from
% ROOT instead; the bank itself is fixed, so this changes no number.
cfg.BankFile = fullfile(ROOT,'data','scenario_bank_WY2024_N10_seed20260920_T8760.mat');
assert(isfile(cfg.BankFile), 'missing the shock bank: %s', cfg.BankFile);

C0 = dataload(ROOT, cfg.Year, cfg.FloorPct, 7);
B  = load(cfg.BankFile,'bank');  bank = B.bank;
C  = C0;  C.zseq = bank.zseq(:,:,job.scen);

% NOTE: this struct is copied verbatim from sweepRunOne.m. If that file
% ever changes, copy it again -- the gate is only meaningful while the two
% expressions agree.
cr = struct('Horizon',cfg.Horizon,'Epsilon',job.eps,'Kappa',job.kappa, ...
    'Theta',cfg.Theta,'MaxHours',cfg.MaxHours,'ScenarioID',job.scen, ...
    'ScenarioSeed',bank.seeds(job.scen),'SolutionMethod',job.method, ...
    'SSHMaxIterations',20,'SSHPhiTolerance',5e-4, ...
    'SSHBoundaryTolerance',1e-5,'SSHLookaheadHours',cfg.Horizon, ...
    'SSHInteriorPoint','bonferroni','ProgressHours',1e9, ...
    'Verbose',false,'ChanceSides',cfg.ChanceSides,'TrackBand',cfg.TrackBand, ...
    'TrackSide',cfg.TrackSide,'PriceMode',cfg.PriceMode, ...
    'SigmaModel',cfg.SigmaModel,'PriceFlatten',cfg.PriceFlatten, ...
    'PriceNoiseSigma',cfg.PriceNoiseSigma,'PriceNoiseRho',cfg.PriceNoiseRho, ...
    'PriceNoiseSeed',cfg.PriceNoiseSeed,'ChanceMode',cfg.ChanceMode, ...
    'WorldKappa',cfg.WorldKappa,'PolicyKappa',cfg.PolicyKappa, ...
    'PowerModel',cfg.PowerModel);

cr.MarginRule = 'chance';          % the identity case -- must change nothing

fprintf('  running realtimeSimulation (about 5 minutes)...\n');
t0 = tic;
out = realtimeSimulation(C, job.policy, cr);
fprintf('  done in %.1f min\n\n', toc(t0)/60);
new = out.summary;

%% ------------------------------------------------------------------ compare
fn = intersect(fieldnames(ref), fieldnames(new));
bad = {};  skipped = {};  nChecked = 0;
fprintf('  %-38s %16s %16s %10s\n','field','stored','rebuilt','rel diff');
for i = 1:numel(fn)
    f = fn{i};
    if ~isempty(regexpi(f, SKIP_PAT, 'once'))
        skipped{end+1} = f; %#ok<SAGROW>
        continue
    end
    a = ref.(f);  b = new.(f);
    if ~(isnumeric(a) || islogical(a)) || ~isscalar(a), continue; end
    if ~(isnumeric(b) || islogical(b)) || ~isscalar(b), continue; end
    a = double(a);  b = double(b);
    nChecked = nChecked + 1;
    if isnan(a) && isnan(b), continue; end
    rel = abs(a-b)/max(abs(a),1);
    if rel > TOL
        bad{end+1} = f; %#ok<SAGROW>
        fprintf('  %-38s %16.6g %16.6g %10.2e  <-- MISMATCH\n', f, a, b, rel);
    end
end

% the added fields should be present and inert
assert(isfield(new,'margin_rule') && strcmp(new.margin_rule,'chance'), ...
    'realtimeSimulation did not record margin_rule = chance');
fprintf('\n  margin_rule = %s, margin_z = %.6f\n', new.margin_rule, new.margin_z);

% Name the excluded fields rather than counting them, so an exclusion that
% should not have happened is visible in the transcript instead of hidden
% behind a number.
fprintf('  %d scalar fields compared exactly\n', nChecked);
fprintf('  %d excluded as timing: %s\n', numel(skipped), strjoin(skipped,', '));
if isempty(bad)
    fprintf('\n  PASS -- realtimeSimulation reproduces the published run at MarginRule = chance.\n');
    fprintf('  The benchmark runs may proceed.\n');
else
    error('verifyParity:mismatch', ...
        ['FAIL -- %d field(s) differ: %s\n' ...
         'The copy is not inert. Diff src/realtimeSimulation.m against ' ...
         'src/realtimeSimulation.m; the patch should be exactly three hunks ' ...
         '(header, margin block after the ChanceMode line, summary field). ' ...
         'Do not run any benchmark until this passes.'], ...
        numel(bad), strjoin(bad,', '));
end

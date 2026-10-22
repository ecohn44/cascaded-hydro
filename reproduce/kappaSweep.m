%% ========================================================================
% REPRODUCE: KAPPA SWEEP
% ========================================================================
%
% Risk attitude of the DDU controller: sweep kappa at the
%                 finalised operating point.
%
%   What kappa does (realtimeSimulation.m ~line 72):
%       gammaModel = cfg.Kappa * gammaData
%       DDU:  gamma = gammaModel
%             sig0  = max(sigDIU - gamma*meanDu, 0.2*sigDIU)
%   so kappa steepens the ramp-sensitive slope while HOLDING the fitted
%   conditional standard deviation at the mean ramp. It is a pure SHAPE
%   parameter: how much of the variance the controller attributes to
%   ramping rather than to baseline. Average conservatism does not move.
%
%   The realised world is untouched (realtimeSimulation.m ~line 611 uses gammaData,
%   never gammaModel), so every kappa faces the identical inflow shocks.
%   This is the controlled comparison -- Test 1 in
%   claude/RESULTS_kappa_scaling_tests.md -- not the level handicap.
%
%   WHY THE GRID STOPS AT 1.50. The sig0 clamp binds when
%       kappa > 0.8 * sigDIU / (gammaData*meanDu)
%   Evaluated on the WY2024 fit:
%       MCN 16.69   JDA 8.28   TDA 1.52   BON 6.22
%   The Dalles is the binding unit at kappa = 1.52, because 52.6% of its
%   sigma is already ramp-attributed (5-13% at the other three). Past that
%   point TDA's sig0 sits on the floor and further kappa RAISES its total
%   sigma instead of redistributing it, so the sweep silently turns into
%   "more conservative at TDA" and stops measuring risk attitude. Do not
%   extend the grid above 1.5 without re-deriving that bound.
%
%   kappa = 0 is not in the grid because runSweep already emits the DIU
%   job at kappa 0 for every (scenario, epsilon), and realtimeSimulation takes the
%   identical arithmetic path there. DIU IS the left-hand anchor.
%
%   Jobs: 4 scenarios x (1 DIU + 6 DDU) = 28, SSH only, ~5.2 min each
%   measured from the headline fill log, so budget about 2.5 hours.
%
%   BEFORE RUNNING: confirm no stale workers are alive --
%       pgrep -fl -i matlab      should come back empty

clear; clc;

%% ---------------------------------------------------------------- CONFIG
THETA = 400;                                   % finalised tracking weight
EPS   = 0.05;                                  % see the note below
KAPPA = [0.25 0.50 0.75 1.00 1.25 1.50];       % -> k025 k050 k075 k100 k125 k150
SCENS = 1:4;

% EPS = 0.05 rather than 0.10. At theta = 400 the DDU margin is wider at
% 0.05 on both axes that this figure plots -- IVI -23.3% vs -18.4%, energy
% +0.041% vs +0.030% -- so the sweep has more room to show a trend, and it
% is the epsilon whose trajectories are already on disk, which keeps this
% figure on the same operating point as the risk-allocation figure.

TAG = sprintf('BPA_KAP_TH%07d_E%03d', THETA, round(1000*EPS));

% This file lives in functions/, so the repository root is two levels up.
logFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    sprintf('logs/kappa_sweep_%s.log', datestr(now,'yyyymmdd_HHMMSS')));
if ~isfolder(fileparts(logFile)), mkdir(fileparts(logFile)); end
diary(logFile); diary on;
fprintf('Log: %s\n', logFile);
fprintf('Started %s\n', datestr(now));
fprintf('tag %s | theta %d | eps %.3f | kappa %s\n', ...
        TAG, THETA, EPS, mat2str(KAPPA));

% Byte-identical to run_ssh_fill.m except that Kappa is now the swept axis.
% WorldKappa stays at 1: the environment must not move with the controller.
base = { ...
    'Horizon',24, 'Year',2024, ...
    'SigmaModel','archx', 'PowerModel','mccormick', ...
    'ChanceMode','on', 'ChanceSides','lower', ...
    'TrackSide','both', 'TrackBand',0, 'FloorPct',0.20, ...
    'PriceMode','bpa', 'PriceFlatten',0, ...
    'PriceNoiseSigma',0.30, 'PriceNoiseRho',0.80, ...
    'WorldKappa',1, 'PolicyKappa',1, ...
    'Workers',7 };

t00 = tic;
for sc = SCENS
    fprintf('\n===== %s : SSH, scenario %d =====\n', TAG, sc);
    t0 = tic;
    try
        runSweep('Tag',TAG, base{:}, ...
            'Theta',THETA, 'Scenarios',sc, 'Epsilon',EPS, ...
            'Kappa',KAPPA, 'Methods',{'SSH'}, 'MaxHours',8760, ...
            'SaveTrajEpsilon',EPS, 'SaveTrajKappa',[0 KAPPA]);
    catch ME
        fprintf('!! scenario %d aborted: %s\n', sc, ME.message);
    end
    fprintf('scenario %d took %.1f min (elapsed %.1f h)\n', ...
            sc, toc(t0)/60, toc(t00)/3600);
end

%% ---------------------------------------------------------------- REPORT
root   = fileparts(fileparts(mfilename('fullpath')));
jobDir = fullfile(root,'results','runs',TAG,'jobs');
expect = numel(SCENS) * (1 + numel(KAPPA));
n      = numel(dir(fullfile(jobDir,'*_SSH.mat')));
fprintf('\n===== kappa sweep =====\n');
fprintf('  SSH jobs %d/%d\n', n, expect);

for kk = [0 KAPPA]
    pol = 'DDU'; if kk == 0, pol = 'DIU'; end
    m = numel(dir(fullfile(jobDir, sprintf('*_%s_k%03d_SSH.mat', pol, round(100*kk)))));
    fprintf('    kappa %4.2f (%s) : %d/%d\n', kk, pol, m, numel(SCENS));
end

errs = dir(fullfile(jobDir,'*_error.txt'));
if isempty(errs)
    if n == expect
        fprintf('\nNo error files. Sweep complete.\n');
    else
        fprintf('\nNo error files, but jobs are missing -- check the log above.\n');
    end
else
    fprintf('\n%d job(s) FAILED:\n', numel(errs));
    for q = 1:numel(errs)
        fid = fopen(fullfile(errs(q).folder, errs(q).name));
        fprintf('  %-38s %s\n', errs(q).name, strtrim(fgetl(fid)));
        fclose(fid);
    end
end

% A clamp warning is the one silent failure mode worth checking by hand:
% if TDA's planned slack collapses at k150 relative to k125, the pivot has
% stopped holding and the top point should be dropped from the figure.
fprintf('\nCheck before plotting: TDA behaviour at k150 vs k125.\n');
fprintf('Total %.1f h. Finished %s\n', toc(t00)/3600, datestr(now));
diary off;

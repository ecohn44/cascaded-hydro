%% ========================================================================
% TEST: SETTLEMENT NOISE IS SCORING-ONLY
% ========================================================================
%
% Two verifications that must hold before a night of solver time is spent.
%
%   Run it:  >> testSettlementNoise
%
%   1. SSH completes under PriceMode 'bpa'. The tracking-weight calibration
%      was BON only, so SSH on the published-tariff price path is otherwise
%      untested.
%   2. Settlement noise does not touch dispatch. The code claims
%      PriceNoiseSigma affects scoring only: every physical metric
%      bit-identical to sigma = 0, and only revenue_usd moving. This proves
%      it rather than assuming it, which matters because the headline grid
%      runs with the noise ON.
%
% The headline grid is only trustworthy if check 2 passes. If it fails, the
% published revenue_usd column is contaminated by a dispatch change and the
% grid must be rerun at PriceNoiseSigma = 0.
%
% RUNTIME. About ten minutes: a 168 h SSH pair and two 720 h BON runs.
%
% Tags BPA_check_ssh, BPA_check_noise0, BPA_check_noise3. These are test
% artefacts, not paper results.

clear; clc;

ROOT = fileparts(fileparts(mfilename('fullpath')));
addpath(ROOT); addpath(fullfile(ROOT,'functions'));

% Shared settings. Methods, Epsilon, Theta, MaxHours, Scenarios, Kappa and
% PriceNoiseSigma are passed per call so nothing is specified twice.
base = { ...
    'Horizon',24, 'Year',2024, ...
    'SigmaModel','archx', 'PowerModel','mccormick', ...
    'ChanceMode','on', 'ChanceSides','lower', ...
    'TrackSide','both', 'TrackBand',0, 'FloorPct',0.20, ...
    'PriceMode','bpa', 'PriceFlatten',0, ...
    'Workers',7 };

%% 1. SSH under the published tariff
fprintf('\n=== Check 1: SSH under PriceMode bpa ===\n');
runSweep('Tag','BPA_check_ssh', base{:}, ...
    'Scenarios',1, 'Kappa',1, 'Methods',{'BON','SSH'}, ...
    'MaxHours',168, 'Epsilon',0.05, 'Theta',400, 'PriceNoiseSigma',0);

%% 2. Settlement noise must not touch dispatch
fprintf('\n=== Check 2: settlement noise is scoring-only ===\n');
runSweep('Tag','BPA_check_noise0', base{:}, ...
    'Scenarios',1, 'Kappa',1, 'Methods',{'BON'}, ...
    'MaxHours',720, 'Epsilon',0.05, 'Theta',400, 'PriceNoiseSigma',0);
runSweep('Tag','BPA_check_noise3', base{:}, ...
    'Scenarios',1, 'Kappa',1, 'Methods',{'BON'}, ...
    'MaxHours',720, 'Epsilon',0.05, 'Theta',400, 'PriceNoiseSigma',0.30);

A = readtable(fullfile(ROOT,'results','runs', ...
    'BPA_check_noise0','summary_BPA_check_noise0.csv'));
B = readtable(fullfile(ROOT,'results','runs', ...
    'BPA_check_noise3','summary_BPA_check_noise3.csv'));

%% 3. Verdict
fprintf('\n--- settlement noise verdict ---\n');
mustMatch = {'energy_MWh','lower_ivi_kaf','step_infeas','spill_kaf', ...
             'lower_exit_pct','revenue_sched_usd'};
ok = true;
for k = 1:numel(mustMatch)
    f = mustMatch{k};
    d = max(abs(A.(f) - B.(f)));
    verdict = 'IDENTICAL';
    if d > 1e-9, verdict = '*** DIFFERS ***'; ok = false; end
    fprintf('  %-20s max |diff| = %.3e   %s\n', f, d, verdict);
end

dRev = max(abs(A.revenue_usd - B.revenue_usd));
if dRev > 1
    revVerdict = 'DIFFERS, as it should';
else
    revVerdict = '*** SHOULD DIFFER ***';
end
fprintf('  %-20s max |diff| = %.3e   %s\n', 'revenue_usd', dRev, revVerdict);

if ok && dRev > 1
    fprintf('\nPASS: dispatch is unchanged by settlement noise; only revenue moves.\n');
else
    error('testSettlementNoise:failed', ...
        ['FAIL: settlement noise changed dispatch. The headline grid runs ' ...
         'with PriceNoiseSigma = 0.30, so its physical metrics cannot be ' ...
         'trusted until this is resolved.']);
end

% price_mean is the mean over the HOURS ACTUALLY RUN, not the annual mean: a
% 720 h window starting 1 Oct covers October, so about 41.1 is correct here.
% Only a full 8760 h run should read 39.64.
fprintf('\nprice_mode must read bpa. price_mean below is over the hours run.\n');
disp(A(:,{'policy','method','price_mode','price_mean','revenue_usd'}));

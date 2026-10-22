%% ========================================================================
% REPRODUCE: HEADLINE GRID
% ========================================================================
%
% The paper's main results grid: DDU vs DIU across epsilon, under both
% Bonferroni and SSH, at three tracking weights.
%
%   Run it:  >> headlineGrid
%
% Produces tags BPA_HL2_TH{0000200,0000400,0000800}. Theta = 400 is the
% calibrated choice and the one the headline table reports; 200 and 800
% bracket it so the choice is shown to be a plateau rather than a knife edge.
%
% Settings are the locked paper configuration, listed verbatim below rather
% than pulled from a shared defaults file, so that what produced the
% published numbers is readable in one place. See README section 6.
%
%   Price   : BPA PF-24 Load Shaping rate, as published, no rescaling.
%             Schedule PF-24 section 2.1.3.1, FY2024-2025, FERC-approved.
%             data/bpa_pf24_lmp_wy2024.csv, annual mean $39.64/MWh.
%   Sigma   : ARCH-X, data/archx_params.csv.
%   Guard   : two-sided water-value guard, one-sided chance constraint, P20.
%   Noise   : PriceNoiseSigma 0.30, rho 0.80. Scoring only -- priceNoiseChecks
%             proves dispatch is bit-identical to sigma = 0 and only
%             revenue_usd moves. revenue_sched_usd is revenue at the
%             published tariff. The seed comes from the scenario seed, so
%             DIU and DDU in a scenario settle against the SAME realised
%             path and stay paired.
%   Kappa   : PolicyKappa = 1, WorldKappa = 1, Kappa = 1, so DDU carries the
%             fitted conditional sigma and the world is the fitted model.
%
% RUNTIME. A full-year BON job measures about 320 s and SSH about 2050 s. At
% 4 scenarios that is 96 jobs per theta, roughly 4.5 h on 7 workers, so the
% three thetas together are an overnight run. runSweep is resumable under the
% same Tag: it skips finished jobs, so a stopped run resumes where it left
% off, and raising Scenarios later re-runs only the new ones.
%
% DO NOT change a physics setting under an existing Tag. Mixed settings in
% one tag folder cannot be told apart afterwards.

clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root,'functions'));

base = { ...
    'Horizon',24, 'Year',2024, ...
    'SigmaModel','archx', 'PowerModel','mccormick', ...
    'ChanceMode','on', 'ChanceSides','lower', ...
    'TrackSide','both', 'TrackBand',0, 'FloorPct',0.20, ...
    'PriceMode','bpa', 'PriceFlatten',0, ...
    'PriceNoiseSigma',0.30, 'PriceNoiseRho',0.80, ...
    'WorldKappa',1, 'Kappa',1, 'PolicyKappa',1, ...
    'Workers',7 };

EPS    = [0.01 0.05 0.10 0.15 0.20 0.25];
SCEN   = 1:4;
THETAS = [200 400 800];        % below / chosen / above

fprintf('headline grid: %d jobs across %d theta values\n', ...
    numel(EPS)*2*2*numel(SCEN)*numel(THETAS), numel(THETAS));
fprintf('  eps    : %s\n', mat2str(EPS));
fprintf('  theta  : %s   (400 is the calibrated choice)\n', mat2str(THETAS));
fprintf('  scen   : %s\n', mat2str(SCEN));
fprintf('  price  : bpa (full year, expect price_mean = 39.64 every row)\n\n');

for th = THETAS
    tag = sprintf('BPA_HL2_TH%07d', th);
    fprintf('\n===== %s : theta = %d =====\n', tag, th);
    runSweep('Tag',tag, base{:}, ...
        'Theta',th, ...
        'Scenarios',SCEN, 'Epsilon',EPS, ...
        'Methods',{'BON','SSH'}, ...
        'MaxHours',8760, ...
        'SaveTrajEpsilon',0.05, 'SaveTrajKappa',[0 1]);
end

fprintf('\nDone. results/runs/BPA_HL2_TH{0000200,0000400,0000800}/\n');
fprintf('Check price_mean is 39.64 on every row before using the results.\n');

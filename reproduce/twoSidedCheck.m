%% ========================================================================
% REPRODUCE: TWO-SIDED CHANCE CONSTRAINT
% ========================================================================
%
% Evidence that restricting the joint chance constraint to the lower storage
% bound is an exact reduction rather than a modelling convenience.
%
%   Run it:  >> twoSidedCheck
%
% Produces tags BPA_2SIDE_P10_TH{0000000,0000100,0000200,0000400}.
%
% THE ARGUMENT. Spill makes the upper storage bound hold pathwise, so the
% 2n-tail and n-tail joint constraints are the same inequality and the
% two-sided Bonferroni is dominated: it spends risk budget on a bound that
% cannot be violated. Splitting epsilon over 2n tails instead of n raises the
% quantile from Phi^-1(1-eps/4) to Phi^-1(1-eps/8) and buys nothing. These
% runs are the empirical side of that; the reduction is derived in the
% manuscript and summarised in docs/realtime.md.
%
% The only setting that differs from the paper's configuration is
% ChanceSides, marked below.
%
% RUNTIME. 4 thetas x 2 policies, BON only, full year: 8 jobs, about 10 min
% on 7 workers.

clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root,'functions'));

base = { ...
    'Horizon',24, 'Year',2024, ...
    'SigmaModel','archx', 'PowerModel','mccormick', ...
    'ChanceMode','on', 'ChanceSides','both', ...   % <-- the only change
    'TrackSide','both', 'TrackBand',0, 'FloorPct',0.10, ...
    'PriceMode','bpa', 'PriceFlatten',0, ...
    'PriceNoiseSigma',0.30, 'PriceNoiseRho',0.80, ...
    'WorldKappa',1, 'Kappa',1, 'PolicyKappa',1, ...
    'Workers',7 };

EPS    = 0.05;
SCEN   = 1;
THETAS = [0 100 200 400];

fprintf('two-sided chance constraint | P10 floor, eps %.2f, BON only\n', EPS);
fprintf('  theta : %s\n', mat2str(THETAS));
fprintf('  z one-sided %.4f  ->  two-sided %.4f\n', ...
        -sqrt(2)*erfcinv(2*(1-EPS/4)), -sqrt(2)*erfcinv(2*(1-EPS/8)));
fprintf('  jobs  : %d\n\n', numel(THETAS)*2);

for th = THETAS
    tag = sprintf('BPA_2SIDE_P10_TH%07d', th);
    fprintf('\n===== %s : theta = %d =====\n', tag, th);
    runSweep('Tag',tag, base{:}, ...
        'Theta',th, ...
        'Scenarios',SCEN, 'Epsilon',EPS, ...
        'Methods',{'BON'}, ...
        'MaxHours',8760, ...
        'SaveTrajEpsilon',0.05, 'SaveTrajKappa',[0 1]);
end

fprintf(['\nDone. Compare against BPA_FLOOR_P10_TH0000400 (one-sided, same\n' ...
         'floor and epsilon) at theta = 400. The number that answers the\n' ...
         'question is BON''s position in its band, not generation.\n']);

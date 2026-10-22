%% ========================================================================
% REPRODUCE: FLOOR PERCENTILE SENSITIVITY
% ========================================================================
%
% Robustness of the DDU result to the choice of lower envelope percentile.
%
%   Run it:  >> floorSensitivity
%
% Produces tags BPA_FLOOR_P10_TH0000400 and BPA_FLOOR_P30_TH0000400. P20 is
% not rerun here: it already exists as the BON rows of BPA_HL2_TH0000400.
%
% READING IT. Compare step_infeas and lower_exit_pct across floors. Do NOT
% compare lower_ivi_kaf across floors: each run measures violation against
% its own envelope, so the three numbers are distances to three different
% targets and the comparison is meaningless.
%
% RUNTIME. 2 floors x 6 eps x 2 policies, BON only, full year: 24 jobs at
% about 320 s, roughly 20 min on 7 workers.

clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root,'functions'));

base = { ...
    'Horizon',24, 'Year',2024, ...
    'SigmaModel','archx', 'PowerModel','mccormick', ...
    'ChanceMode','on', 'ChanceSides','lower', ...
    'TrackSide','both', 'TrackBand',0, ...
    'PriceMode','bpa', 'PriceFlatten',0, ...
    'PriceNoiseSigma',0.30, 'PriceNoiseRho',0.80, ...
    'WorldKappa',1, 'Kappa',1, 'PolicyKappa',1, ...
    'Workers',7 };

EPS    = [0.01 0.05 0.10 0.15 0.20 0.25];
SCEN   = 1;                    % one scenario: directional, no intervals
THETA  = 400;
FLOORS = [0.10 0.30];          % P20 already exists as BPA_HL2_TH0000400

fprintf('floor sensitivity at theta = %d, BON only\n', THETA);
fprintf('  floors : %s   (P20 is the existing BPA_HL2_TH0000400)\n', mat2str(FLOORS));
fprintf('  eps    : %s\n', mat2str(EPS));
fprintf('  scen   : %s\n', mat2str(SCEN));
fprintf('  jobs   : %d\n\n', numel(FLOORS)*numel(EPS)*2*numel(SCEN));

for fp = FLOORS
    tag = sprintf('BPA_FLOOR_P%02d_TH%07d', round(100*fp), THETA);
    fprintf('\n===== %s : floor P%02.0f =====\n', tag, 100*fp);
    runSweep('Tag',tag, base{:}, ...
        'Theta',THETA, 'FloorPct',fp, ...
        'Scenarios',SCEN, 'Epsilon',EPS, ...
        'Methods',{'BON'}, ...
        'MaxHours',8760, ...
        'SaveTrajEpsilon',0.05, 'SaveTrajKappa',[0 1]);
end

fprintf(['\nDone. Compare against BPA_HL2_TH0000400 (P20, BON rows only).\n' ...
         'Read step_infeas and lower_exit_pct across floors; do NOT compare\n' ...
         'lower_ivi_kaf across floors, since each is measured against its own\n' ...
         'envelope.\n']);

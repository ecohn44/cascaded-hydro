%% ========================================================================
% REPRODUCE: TRACKING-WEIGHT SELECTION
% ========================================================================
%
% The theta selection panel: how the tracking weight trades violation
% against generation, and why theta = 400 is the choice.
%
%   Run it:  >> thetaSelection
%
% Produces tags BPA_TH7_P20_TH* at ten paired scenarios, BON only.
%
% WHAT THIS PANEL RESOLVES, AND WHAT IT DOES NOT. The selection resolves
% cleanly on hours below target and not at all on IVI:
%
%     hours   residual sd  13.4 h     F(6,30) = 53.35   p = 1.1e-14
%     IVI     residual sd  90.6 kaf   F(6,30) =  0.89   p = 0.51
%
% Hours counts events while IVI integrates depth, so one harsh inflow draw
% moves IVI far more than it moves the count. That is a property of the
% metric, not a defective scenario: guard, fallback and emergency counts are
% all zero across scenarios, and spill and step-infeasibility are flat. Ten
% scenarios is the whole shock bank, drawn from one seed so every theta is
% compared on the SAME inflow paths. Report hours as the selection evidence
% and treat IVI at a single theta as indicative only.
%
% RUNTIME. 7 theta values x 10 scenarios x 2 policies, BON only, full year:
% 140 jobs at about 320 s, roughly 1.8 h on 7 workers.

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

SPAN = [200 400 600 800 1000 1200 1400];
EPS  = 0.05;
SCEN = 1:10;

fprintf('theta selection, P20, eps %.2f, BON only\n', EPS);
fprintf('  theta     : %s\n', mat2str(SPAN));
fprintf('  scenarios : %s   (the whole shock bank, one seed)\n', mat2str(SCEN));
fprintf('  jobs      : %d\n\n', numel(SPAN)*2*numel(SCEN));

for th = SPAN
    tag = sprintf('BPA_TH7_P20_TH%07d', th);
    fprintf('\n===== %s =====\n', tag);
    runSweep('Tag',tag, base{:}, ...
        'Theta',th, 'Scenarios',SCEN, 'Epsilon',EPS, ...
        'Methods',{'BON'}, 'MaxHours',8760, ...
        'SaveTrajEpsilon',0.05, 'SaveTrajKappa',[0 1]);
end

fprintf('\nDone. Tags BPA_TH7_P20_*.\n');
fprintf('Render with figures/make_fig_theta_select_hours.py.\n');

%% ========================================================================
% PAPER RESULTS MANIFEST
% ========================================================================

function M = paperTables()
% The single source of truth linking each paper result to the code behind it
%
% Outputs:
%   M  - struct array, one element per published table or figure, with
%        fields: name, tags, driver, render, artefact, jobs, hours, reproduce
%
% Every published number in the paper traces through one row of this
% manifest: the tag whose summary CSV holds it, the driver script that
% produced the tag, and the script that renders it. reproduceTables.m reruns
% the rows marked reproduce = true. Keep it in step with the reproduction
% table in README section 6, which lists the same results.
%
% FIELDS.
%   name       - what the result is called in the paper
%   tags       - the result tags it reads, as a cell array
%   driver     - the script that produces those tags
%   render     - the script that turns them into the table or figure
%   artefact   - the shipped file whose presence proves the result is here,
%                relative to the repository root. Empty means the result is
%                reproduced on demand and nothing is shipped for it.
%   jobs       - number of solver jobs behind it
%   hours      - rough wall-clock on 7 workers
%   reproduce  - true if reproduceTables.m reruns it by default
%
% reproduce = false does not mean unverified. It means the row is too
% expensive to rerun on every pass and is covered instead by the parity gate
% and the shipped CSV. Set it true to include the row.

k = 0;

k = k + 1;
M(k).name      = 'Headline grid: DDU vs DIU over epsilon, BON and SSH';
M(k).tags      = {'BPA_HL2_TH0000400'};
M(k).driver    = 'reproduce/headlineGrid.m';
M(k).render    = 'reproduce/makeTables.m, figures/make_ddu_heatmaps.m, figures/make_ssh_vs_bon.m';
M(k).artefact  = 'results/runs/BPA_HL2_TH0000400/summary_BPA_HL2_TH0000400.csv';
M(k).jobs      = 96;
M(k).hours     = 4.5;
M(k).reproduce = true;

k = k + 1;
M(k).name      = 'Decision-dependence sweep over kappa';
M(k).tags      = {'BPA_KAP_TH0000400_E050'};
M(k).driver    = 'reproduce/kappaSweep.m';
M(k).render    = 'figures/make_kappa_figure.m';
M(k).artefact  = 'results/runs/BPA_KAP_TH0000400_E050/summary_BPA_KAP_TH0000400_E050.csv';
M(k).jobs      = 28;
M(k).hours     = 1.6;
M(k).reproduce = true;

k = k + 1;
M(k).name      = 'Tracking-weight stability at theta 100, 200 and 800';
% make_ddu_heatmaps and make_tracking_figure both read theta = 100, so it is
% part of the published ladder even though 400 is the headline. Listing it here
% is what makes checkEnvironment notice if it ever goes missing again.
M(k).tags      = {'BPA_HL2_TH0000100', 'BPA_HL2_TH0000200', 'BPA_HL2_TH0000800'};
M(k).driver    = 'reproduce/headlineGrid.m';
M(k).render    = 'figures/make_tracking_figure.m';
M(k).artefact  = 'results/runs/BPA_HL2_TH0000800/summary_BPA_HL2_TH0000800.csv';
M(k).jobs      = 288;
M(k).hours     = 13.5;
M(k).reproduce = false;

k = k + 1;
M(k).name      = 'Tracking-weight selection over theta';
M(k).tags      = {'BPA_TH7_P20_TH0000200','BPA_TH7_P20_TH0000400', ...
                  'BPA_TH7_P20_TH0000600','BPA_TH7_P20_TH0000800', ...
                  'BPA_TH7_P20_TH0001000','BPA_TH7_P20_TH0001200', ...
                  'BPA_TH7_P20_TH0001400'};
M(k).driver    = 'reproduce/thetaSelection.m';
M(k).render    = 'figures/make_fig_theta_select_hours.py';
M(k).artefact  = 'results/runs/BPA_TH7_P20_TH0000400/summary_BPA_TH7_P20_TH0000400.csv';
M(k).jobs      = 140;
M(k).hours     = 1.8;
M(k).reproduce = false;

k = k + 1;
M(k).name      = 'Deterministic and robust MPC benchmarks';
M(k).tags      = {};              % written straight to results/runs/benchmarks.csv
M(k).driver    = 'reproduce/benchmarks.m';
M(k).render    = 'reproduce/makeTables.m, reproduce/analyzeBenchmarks.m';
M(k).artefact  = 'results/runs/benchmarks.csv';
M(k).jobs      = 28;
M(k).hours     = 1.5;
M(k).reproduce = false;

k = k + 1;
M(k).name      = 'Floor percentile sensitivity, P10 and P30';
M(k).tags      = {'BPA_FLOOR_P10_TH0000400', 'BPA_FLOOR_P30_TH0000400'};
M(k).driver    = 'reproduce/floorSensitivity.m';
M(k).render    = 'none: read from the summary CSVs';
M(k).artefact  = 'results/runs/BPA_FLOOR_P30_TH0000400/summary_BPA_FLOOR_P30_TH0000400.csv';
M(k).jobs      = 24;
M(k).hours     = 0.35;
M(k).reproduce = false;

k = k + 1;
M(k).name      = 'Two-sided chance constraint comparison';
M(k).tags      = {'BPA_2SIDE_P10_TH0000000','BPA_2SIDE_P10_TH0000100', ...
                  'BPA_2SIDE_P10_TH0000200','BPA_2SIDE_P10_TH0000400'};
M(k).driver    = 'reproduce/twoSidedCheck.m';
M(k).render    = 'none: read from the summary CSVs';
M(k).artefact  = 'results/runs/BPA_2SIDE_P10_TH0000400/summary_BPA_2SIDE_P10_TH0000400.csv';
M(k).jobs      = 8;
M(k).hours     = 0.2;
M(k).reproduce = false;

k = k + 1;
M(k).name      = 'Scalability timing at 4, 8 and 12 units';
M(k).tags      = {};              % timing only; nothing is shipped for it
M(k).driver    = 'reproduce/scalability.m';
M(k).render    = 'figures/make_ssh_vs_bon.m';
M(k).artefact  = '';
M(k).jobs      = 12;
M(k).hours     = 2.0;
M(k).reproduce = false;

end

%MAKE_DDU_HEATMAPS  DDU versus DIU across the theta x epsilon grid, SSH.
%
%   Run it:  >> make_ddu_heatmaps        (script, not a function)
%
%   Replaces the heatmap section of baseline_figures.m, which was built on a
%   kappa x epsilon grid read from paper/plot_data/14_tracking_toggle.csv.
%   The final sweep holds kappa fixed (DIU = 0, DDU = 1) and varies theta, so
%   the axes are theta x epsilon and the data is read straight from the four
%   BPA_HL2 summary files -- no intermediate CSV to go stale.
%
%   Every cell is a paired DDU - DIU comparison, matched by scenario, then
%   averaged over the four scenarios. Both metrics are signed so that
%   POSITIVE ALWAYS MEANS DDU IS BETTER:
%       generation  (DDU - DIU)/1e3, GWh     (more energy)
%       violation   100*(DIU - DDU)/DIU      (less IVI)
%   The colour scale spans the observed range, min to max, so the full ramp
%   is used to separate cells rather than compressed against one end. The
%   metric name labels the colorbar instead of sitting above the axes.
%
%   Outputs, two separate figures:
%       figures_out/ddu_heatmap_generation.{pdf,png}
%       figures_out/ddu_heatmap_violation.{pdf,png}

clear; clc; close all;

%% ---------------------------------------------------------------- CONFIG
here = fileparts(mfilename('fullpath'));
ROOT = fileparts(here);
FIGD = fullfile(ROOT,'figures_out');
if ~isfolder(FIGD), mkdir(FIGD); end

THETAS  = [800 400 200 100];                     % largest at top
EPSVALS = [0.25 0.20 0.15 0.10 0.05 0.01];       % confidence 75% -> 99%
METHOD  = "SSH";
FS      = 14;
PAD     = 0.15;     % extend the scale this far below the smallest cell, as a
                    % fraction of the observed range, so no cell renders white

confLabels = "" + (100*(1-EPSVALS)) + "%";

%% ------------------------------------------------------------ LOAD DATA
T = table();
for th = THETAS
    f = fullfile(ROOT,'results','runs', ...
                 sprintf('BPA_HL2_TH%07d',th), ...
                 sprintf('summary_BPA_HL2_TH%07d.csv',th));
    assert(isfile(f), 'missing summary: %s', f);
    T = [T; readtable(f)]; %#ok<AGROW>
end
T.policy = string(T.policy);
T.method = string(T.method);
T = T(T.method == METHOD, :);
fprintf('%s: %d rows, theta %s\n', METHOD, height(T), mat2str(unique(T.theta)'));

%% ------------------------------------------------- paired DDU - DIU cells
nR = numel(THETAS); nC = numel(EPSVALS);
genChange = nan(nR,nC);
iviChange = nan(nR,nC);
cellN     = nan(nR,nC);

for r = 1:nR
    for c = 1:nC
        sel = abs(T.theta - THETAS(r)) < 1e-8 & ...
              abs(T.eps   - EPSVALS(c)) < 1e-8;
        diu = T(sel & T.policy == "DIU", :);
        ddu = T(sel & T.policy == "DDU", :);

        % match DDU and DIU observations by scenario
        [~, iDDU, iDIU] = intersect(ddu.scen, diu.scen);
        ddu = ddu(iDDU,:);  diu = diu(iDIU,:);
        cellN(r,c) = height(ddu);
        if isempty(ddu), continue; end

        genChange(r,c) = mean((ddu.energy_MWh - diu.energy_MWh)/1e3);
        iviChange(r,c) = mean(100*(diu.lower_ivi_kaf - ddu.lower_ivi_kaf) ...
                                 ./ diu.lower_ivi_kaf);
    end
end

if any(cellN(:) ~= 4)
    warning('not every cell has n=4; counts range %d to %d', ...
            min(cellN(:)), max(cellN(:)));
end
if any(genChange(:) < 0) || any(iviChange(:) < 0)
    warning(['at least one cell favours DIU; with a min-to-max scale the ' ...
             'lightest cell is no longer necessarily zero -- read the ' ...
             'colorbar, not the shade, and check the printed ranges below']);
end

%% ---------------------------------------------------------------- FIGURES
% sequential single hue, white at zero through to the largest value
BLUE = [0 114 178]/255;
g    = linspace(0,1,256)';
cmap = 1 - g.*(1 - BLUE);

panels = {genChange, iviChange};
pname  = {'Generation gain (GWh)', 'Violation reduction (-%)'};
pfile  = {'ddu_heatmap_generation', 'ddu_heatmap_violation'};
pfmt   = {'%.3f', '%.1f'};

for p = 1:2
    A  = panels{p};
    mn = min(A(:),[],'omitnan');
    hi = max(A(:),[],'omitnan');
    lo = mn - PAD*(hi - mn);           % floor sits below the smallest cell,
                                       % so even that one carries visible tint

    fig = figure('Color','w','Units','inches','Position',[1 1 6.6 3.7]);
    ax  = axes(fig);
    imagesc(ax, A, [lo hi]); colormap(ax, cmap); hold(ax,'on');

    for r = 1:nR
        for c = 1:nC
            v = A(r,c); if isnan(v), continue; end
            shade = (v - lo)/(hi - lo);
            txt = [0 0 0] + (shade > 0.55)*([1 1 1]-[0 0 0]);
            text(ax, c, r, sprintf(pfmt{p}, v), 'HorizontalAlignment','center', ...
                 'FontSize', FS-1, 'FontName','Times New Roman', 'Color', txt);
        end
    end

    xticks(ax,1:nC); xticklabels(ax, confLabels);
    yticks(ax,1:nR); yticklabels(ax, "\theta = " + string(THETAS));
    xlabel(ax,'Chance-constraint confidence 1 - \epsilon');
    ax.FontName = 'Times New Roman'; ax.FontSize = FS;
    ax.TickLength = [0 0];
    cb = colorbar(ax); cb.FontSize = FS-1;
    cb.Label.String = pname{p};        % the panel title now labels the bar
    cb.Label.FontSize = FS+1;

    exportgraphics(fig, fullfile(FIGD,[pfile{p} '.pdf']), 'ContentType','vector');
    exportgraphics(fig, fullfile(FIGD,[pfile{p} '.png']), 'Resolution', 300);
end

%% ---------------------------------------------------------------- REPORT
fprintf('  generation gain     %+.3f to %+.3f GWh (median %+.3f)\n', ...
    min(genChange,[],'all'), max(genChange,[],'all'), median(genChange,'all'));
fprintf('  violation reduction %+.1f to %+.1f %%   (median %+.1f)\n', ...
    min(iviChange,[],'all'), max(iviChange,[],'all'), median(iviChange,'all'));
fprintf('  cells favouring DDU: generation %d/%d, violation %d/%d\n', ...
    sum(genChange(:) > 0), numel(genChange), ...
    sum(iviChange(:) > 0), numel(iviChange));
fprintf('wrote %s and %s to %s\n', pfile{1}, pfile{2}, FIGD);

%% Plot Driver for DIU vs DDU Comparison (All Units)
% ========================================================================
% Author: Eliza Cohn
% Description: Overlays policy behavior derived during optimization sims
%              for all units. For each unit k, loads that unit's DET, DIU,
%              and DDU optimization results separately and produces:
%                 1) Inflow std dev (M2 vs M3)
%                 2) SoC trajectory + effective bounds (M1, M2, M3)
%                 3) Power dispatch (M1, M2, M3)
%              At the end, aggregates energy across all units.
% ========================================================================

clc; clear; close all;

resultsPath = './resultsBonferroni/';
if ~exist(resultsPath,'dir')
    error('results folder not found: %s', resultsPath);
end

% Figure / font settings
fontAxes   = 14;
fontLegend = 13;
fontTitle  = 16;

% ---- Colors (define once so we can reuse for the final figure) ----
blueBound  = [0 0 1];
greenBound = [0 1 0];
redBound   = [1 0 0];

detFiles = dir(fullfile(resultsPath, 'results_unit*_det.mat'));
if isempty(detFiles)
    error('No results_unit*_det.mat files found in %s', resultsPath);
end

% Parse unit indices from filenames
unitIdx = zeros(numel(detFiles),1);
for k = 1:numel(detFiles)
    fname = detFiles(k).name;               % e.g. 'results_unit3_det.mat'
    tokens = regexp(fname, 'results_unit(\d+)_det\.mat', 'tokens', 'once');
    unitIdx(k) = str2double(tokens{1});
end
unitIdx = sort(unitIdx);
n_units = numel(unitIdx);

fprintf('Found %d units: %s\n', n_units, mat2str(unitIdx));

% -------- NEW: accumulators for total power over all units -----------
P_M1_allUnits = [];   % will size after first unit
P_M2_allUnits = [];
P_M3_allUnits = [];

% Loop over each unit
for uu = unitIdx(:)'
    
    fprintf('\n Processing Unit %d \n', uu);

    % Load optimal trajectories for this unit (separate files)
    detFile = fullfile(resultsPath, sprintf('results_unit%d_det.mat', uu));
    diuFile = fullfile(resultsPath, sprintf('results_unit%d_diu.mat', uu));
    dduFile = fullfile(resultsPath, sprintf('results_unit%d_ddu.mat', uu));

    M1 = load(detFile);   % Deterministic
    M2 = load(diuFile);   % DIU
    M3 = load(dduFile);   % DDU

    % Basic time/index info
    season = M1.season;
    s      = M1.sysparams;
    T      = M1.T;
    tt     = (1:T)';

    % If first unit, initialize the total-power arrays
    if isempty(P_M1_allUnits)
        P_M1_allUnits = zeros(T,1);
        P_M2_allUnits = zeros(T,1);
        P_M3_allUnits = zeros(T,1);
    end

    % Indices for this unit in X and V_eff
    % X columns: [V, p, u, s, q] per unit
    baseX   = (uu-1)*5;           % offset for this unit
    colV    = baseX + 1;          % volume column
    colP    = baseX + 2;          % power column

    % V_eff columns: [V1_max V1_min V2_max V2_min ...]
    colVmax = 2*uu - 1;
    colVmin = 2*uu;

    % --------------------------------------------------------------------
    % FIGURE 1: Standard Deviation Comparison for M2 and M3
    % --------------------------------------------------------------------
    f1 = figure('Position',[100 100 850 350]); hold on; grid on;
    plot(tt, M2.std_hat(:,uu), 'LineWidth',3, 'DisplayName','M2 \sigma_t');
    plot(tt, M3.std_hat(:,uu), 'LineWidth',3, 'DisplayName','M3 \sigma_t');
    ylabel('Standard Deviation');
    xlabel('Time (h)')
    title(sprintf('Unit %02d Inflow Uncertainty (%s)', uu, season));
    legend('Location','northeast');
    set(gca, 'FontSize', fontAxes);
    set(findall(gcf,'Type','text'), 'FontSize', fontAxes);

    saveas(f1, fullfile(resultsPath, sprintf('u%d_sigma_diu_vs_ddu.png', uu)));

    % --------------------------------------------------------------------
    % FIGURE 2: Effective Volume Boundaries + Trajectories
    % --------------------------------------------------------------------
    % Extract bounds & trajectories for this unit
    Vmax_M1 = M1.V_eff(:,colVmax);
    Vmin_M1 = M1.V_eff(:,colVmin);
    V_M1    = M1.X(:,colV);

    Vmax_M2 = M2.V_eff(:,colVmax);
    Vmin_M2 = M2.V_eff(:,colVmin);
    V_M2    = M2.X(:,colV);

    Vmax_M3 = M3.V_eff(:,colVmax);
    Vmin_M3 = M3.V_eff(:,colVmin);
    V_M3    = M3.X(:,colV);

    blueTraj   = 0.6 * blueBound;
    greenTraj  = 0.6 * greenBound;
    redTraj    = 0.6 * redBound;

    f2 = figure('Position',[100 480 1100 480]); hold on; grid on;

    % M1: Deterministic bounds and policy
    plot(tt, Vmin_M1, 'LineWidth',3, 'LineStyle','-.', 'Color',blueBound,  ...
        'DisplayName','M1^{Bounds}');
    plot(tt, Vmax_M1, 'LineWidth',3, 'LineStyle','-.', 'Color',blueBound,  ...
        'HandleVisibility','off');
    plot(tt, V_M1,    'LineWidth',2, 'Color',blueTraj, ...
        'DisplayName','M1^{Trajectory}');

    % M2: DIU bounds and policy
    plot(tt, Vmin_M2, 'LineWidth',3, 'LineStyle','-.', 'Color',greenBound, ...
        'DisplayName','M2^{Bounds}');
    plot(tt, Vmax_M2, 'LineWidth',3, 'LineStyle','-.', 'Color',greenBound, ...
        'HandleVisibility','off');
    plot(tt, V_M2,    'LineWidth',2, 'Color',greenTraj, ...
        'DisplayName','M2^{Trajectory}');

    % M3: DDU bounds and policy
    plot(tt, Vmin_M3, 'LineWidth',3, 'LineStyle','-.', 'Color',redBound,   ...
        'DisplayName','M3^{Bounds}');
    plot(tt, Vmax_M3, 'LineWidth',3, 'LineStyle','-.', 'Color',redBound,   ...
        'HandleVisibility','off');
    plot(tt, V_M3,    'LineWidth',2, 'Color',redTraj, ...
        'DisplayName','M3^{Trajectory}');

    ylabel('SoC Bounds');
    xlabel('Time (h)')
    set(gca, 'FontSize', fontAxes);
    set(findall(gcf,'Type','text'), 'FontSize', fontAxes);
    legend('Location','best', 'FontSize', fontLegend);
    title(sprintf('Unit %02d SoC Trajectory and Effective Bounds', uu), ...
          'FontSize', fontTitle);
    ylim([0 0.5]);

    saveas(f2, fullfile(resultsPath, sprintf('u%d_soc_diu_vs_ddu.png', uu)));

    % --------------------------------------------------------------------
    % FIGURE 3: Power Production (Normalized) - per unit
    % --------------------------------------------------------------------
    p_M1 = M1.X(:,colP);
    p_M2 = M2.X(:,colP);
    p_M3 = M3.X(:,colP);

    % ---- accumulate total power across units for each method ----
    P_M1_allUnits = P_M1_allUnits + p_M1;
    P_M2_allUnits = P_M2_allUnits + p_M2;
    P_M3_allUnits = P_M3_allUnits + p_M3;

    P_all = [p_M1, p_M2, p_M3];

    f3 = figure('Position',[100 480 1100 480]); hold on; grid on;

    b = bar(tt, P_all, 'grouped');

    b(1).FaceColor = blueBound;
    b(2).FaceColor = greenBound;
    b(3).FaceColor = redBound;

    ylabel('Normalized Power (p.u.)');
    xlabel('Time (h)')
    set(gca, 'FontSize', fontAxes);
    set(findall(gcf,'Type','text'), 'FontSize', fontAxes);
    legend({'M1', 'M2', 'M3'}, 'Location', 'best', 'FontSize', fontLegend);
    title(sprintf('Unit %02d Normalized Power Dispatch', uu), ...
          'FontSize', fontTitle);

    saveas(f3, fullfile(resultsPath, sprintf('u%d_power_diu_vs_ddu.png', uu)));

    % (Optional) close figures to avoid clutter
    % close([f1 f2 f3]);
end

% =======================================================================
% Figure 4: Accumulated Energy Over Time 
% =======================================================================

dt = 1;  % hours per time step (adjust if different)
T_all = length(P_M1_allUnits);
tt_all = (1:T_all)';

E_M1 = cumsum(P_M1_allUnits) * dt;
E_M2 = cumsum(P_M2_allUnits) * dt;
E_M3 = cumsum(P_M3_allUnits) * dt;

E_all = [E_M1, E_M2, E_M3];

f4 = figure('Position',[150 300 1100 480]); hold on; grid on;

bE = bar(tt_all, E_all, 'grouped');
bE(1).FaceColor = blueBound;
bE(2).FaceColor = greenBound;
bE(3).FaceColor = redBound;

ylabel('Accumulated Energy (p.u.·h)');
xlabel('Time (h)');
set(gca, 'FontSize', fontAxes);
set(findall(gcf,'Type','text'), 'FontSize', fontAxes);
legend({'M1', 'M2', 'M3'}, 'Location', 'best', 'FontSize', fontLegend);
title('Accumulated Energy Dispatch', 'FontSize', fontTitle);

saveas(f4, fullfile(resultsPath, 'allUnits_accumulated_energy.png'));

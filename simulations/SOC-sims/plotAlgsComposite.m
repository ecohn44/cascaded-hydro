%% Plot Driver for SSH vs BON Comparison 
% ========================================================================
% Author: Eliza Cohn
% Description: Overlays policy behavior derived under BON and SSH algs
%              - One combined 4x2 figure:
%                  * Left column: SoC bounds + trajectories (M1, M2, M3)
%                  * Right column: inflow std dev (M1=0, M2, M3)
%              - One combined accumulated-energy figure (all units).
%
% Assumptions:
%   - X layout per unit: [V, p, u, s, q]
%   - V_eff layout: [V1_max V1_min V2_max V2_min ...]
% ========================================================================

% clc; clear; close all;

%% Plot Settings

bounds_plot = "soc"; %options: "soc" "head"

resultsPath = './resultsBonferroni/wet/';
if ~exist(resultsPath,'dir')
    error('results folder not found: %s', resultsPath);
end

% Figure and font settings
fontAxes   = 14;
fontLegend = 13;
fontTitle  = 16;

% Base colors (used consistently across all figures)
blueBound  = [0 0 1];
greenBound = [0 1 0];
redBound   = [1 0 0];
blackBound = [0 0 0];

% Load files 
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

% Aggregate across units
P_M1_allUnits = [];   % total power over all units (time x 1)
P_M2_allUnits = [];
P_M3_allUnits = [];

% Combined std-dev and bounds figures
std_M2_all  = [];
std_M3_all  = [];
Vmin_M1_all = [];
Vmax_M1_all = [];
V_M1_all    = [];
Vmin_M2_all = [];
Vmax_M2_all = [];
V_M2_all    = [];
Vmin_M3_all = [];
Vmax_M3_all = [];
V_M3_all    = [];


%% Collect Trajectory + Simulation Data
unitCounter = 0;

for uu = unitIdx(:)'   % uu is the actual unit number in the files
    
    unitCounter = unitCounter + 1;   % column index for combined arrays
    idx = unitCounter;
    
    fprintf('\n Processing Unit %d \n', uu);

    % Load optimal trajectories for this unit (separate files)
    detFile = fullfile(resultsPath, sprintf('results_unit%d_det.mat', uu));
    diuFile = fullfile(resultsPath, sprintf('results_unit%d_diu.mat', uu));
    dduFile = fullfile(resultsPath, sprintf('results_unit%d_ddu.mat', uu));

    M1 = load(detFile);   % Deterministic
    M2 = load(diuFile);   % DIU
    M3 = load(dduFile);   % DDU

    % Basic time/index info
    T      = M1.T;
    tt     = (1:T)';

    % Initialize aggregation arrays the first time through
    if isempty(P_M1_allUnits)
        P_M1_allUnits = zeros(T,1);
        P_M2_allUnits = zeros(T,1);
        P_M3_allUnits = zeros(T,1);
        
        % Preallocate combined matrices (T x n_units)
        std_M2_all  = zeros(T, n_units);
        std_M3_all  = zeros(T, n_units);
        Vmin_M1_all = zeros(T, n_units);
        Vmax_M1_all = zeros(T, n_units);
        V_M1_all    = zeros(T, n_units);
        Vmin_M2_all = zeros(T, n_units);
        Vmax_M2_all = zeros(T, n_units);
        V_M2_all    = zeros(T, n_units);
        Vmin_M3_all = zeros(T, n_units);
        Vmax_M3_all = zeros(T, n_units);
        V_M3_all    = zeros(T, n_units);
    end

    % X columns: [V, p, u, s, q] per unit
    baseX   = (uu-1)*5;           % offset for this unit
    colV    = baseX + 1;          % volume column
    colP    = baseX + 2;          % power column

    % V_eff columns indices 
    colVmax = 2*uu - 1;
    colVmin = 2*uu;

    % Standard Deviation 
    std_M2_all(:,idx) = M2.std_hat(:,uu);
    std_M3_all(:,idx) = M3.std_hat(:,uu);

    % SOC Bounds
    Vmax_M1 = M1.V_eff(:,colVmax);
    Vmin_M1 = M1.V_eff(:,colVmin);
    V_M1    = M1.X(:,colV);

    Vmax_M2 = M2.V_eff(:,colVmax);
    Vmin_M2 = M2.V_eff(:,colVmin);
    V_M2    = M2.X(:,colV);

    Vmax_M3 = M3.V_eff(:,colVmax);
    Vmin_M3 = M3.V_eff(:,colVmin);
    V_M3    = M3.X(:,colV);

    Vmin_M1_all(:,idx) = Vmin_M1;
    Vmax_M1_all(:,idx) = Vmax_M1;
    V_M1_all(:,idx)    = V_M1;

    Vmin_M2_all(:,idx) = Vmin_M2;
    Vmax_M2_all(:,idx) = Vmax_M2;
    V_M2_all(:,idx)    = V_M2;

    Vmin_M3_all(:,idx) = Vmin_M3;
    Vmax_M3_all(:,idx) = Vmax_M3;
    V_M3_all(:,idx)    = V_M3;

    % Power
    p_M1 = M1.X(:,colP);
    p_M2 = M2.X(:,colP);
    p_M3 = M3.X(:,colP);

    P_M1_allUnits = P_M1_allUnits + p_M1;
    P_M2_allUnits = P_M2_allUnits + p_M2;
    P_M3_allUnits = P_M3_allUnits + p_M3;

    % Hydraulic Head Coefficients 
    a_M1 = M1.sysparams.a;
    b_M1 = M1.sysparams.b;
    a_M2 = M2.sysparams.a;
    b_M2 = M2.sysparams.b;
    a_M3 = M3.sysparams.a;
    b_M3 = M3.sysparams.b;


end

%% Figure 1: SOC Bounds + Standard Deviation
% Combined "super" figure for bounds + std (4 x 2)
T_all = size(P_M1_allUnits,1);
tt_all = (1:T_all)';

fCombined = figure('Position',[200 50 1400 900]);
tiledlayout(n_units, 2, 'TileSpacing','compact', 'Padding','compact');

% To capture legend handles
hM1_traj = [];
hM2_traj = [];
hM3_traj = [];
hM1_std = [];
hM2_std  = [];
hM3_std  = [];


for i = 1:n_units

    if bounds_plot == "soc"
        V1min = Vmin_M1_all(:,i);
        V1max = Vmax_M1_all(:,i);
        V2min = Vmin_M2_all(:,i);
        V2max = Vmax_M2_all(:,i);
        V3min = Vmin_M3_all(:,i);
        V3max = Vmax_M3_all(:,i);
        ylim_array = [0.9 1.05];%[-.01 0.045];
        vol_title = 'Unit %d SoC Lower Bounds';
        vol_ylabel = 'SoC';
    else
        V1min = a_M1*(Vmin_M1_all(:,i).^b_M1);
        V1max = a_M1*(Vmax_M1_all(:,i).^b_M1);
        V2min = a_M2*(Vmin_M2_all(:,i).^b_M2);
        V2max = a_M2*(Vmax_M2_all(:,i).^b_M2);
        V3min = a_M3*(Vmin_M3_all(:,i).^b_M3);
        V3max = a_M3*(Vmax_M3_all(:,i).^b_M3);
        ylim_array = [-.01 5];
        vol_title = 'Unit %d Hydraulic Head Bounds';
        vol_ylabel = 'Head Height';
    end 

    % Fig 1A (LEFT): SoC bounds + trajectories 
    nexttile; hold on; grid on;
    
    % M1
    p1a = plot(tt_all, V1min, 'LineWidth',3, 'LineStyle','-.', 'Color', blackBound);
    p1b = plot(tt_all, V1max, 'LineWidth',3, 'LineStyle','-.', 'Color', blackBound, 'HandleVisibility','off');
    % h1  = plot(tt_all, V_M1_all(:,i),    'LineWidth',2, 'Color', 0.6*blueBound);
    
    % M2
    p2a = plot(tt_all, V2min, 'LineWidth',3, 'LineStyle','-.', 'Color', blueBound);
    p2b = plot(tt_all, V2max, 'LineWidth',3, 'LineStyle','-.', 'Color', blueBound, 'HandleVisibility','off');
    % h2  = plot(tt_all, V_M2_all(:,i),    'LineWidth',2, 'Color', 0.6*blueBound);
    
    % M3
    p3a = plot(tt_all, V3min, 'LineWidth',3, 'LineStyle','-.', 'Color', redBound);
    p3b = plot(tt_all, V3max, 'LineWidth',3, 'LineStyle','-.', 'Color', redBound, 'HandleVisibility','off');
    % h3  = plot(tt_all, V_M3_all(:,i),    'LineWidth',2, 'Color', 0.6*redBound);
    
    if i == 1
        hM1_traj = p1a; %h1;
        hM2_traj = p2a; %h2;
        hM3_traj = p3a; %h3;
    else
        %set(h1,'HandleVisibility','off');
        %set(h2,'HandleVisibility','off');
        %set(h3,'HandleVisibility','off');
        set(p1a,'HandleVisibility','off');
        set(p2a,'HandleVisibility','off');
        set(p3a,'HandleVisibility','off');
    end
 
    
    xlim([1, T])
    ylim(ylim_array);
    title(sprintf(vol_title, unitIdx(i)), 'FontSize', fontTitle);
    xlabel('Time (h)');  
    ylabel(vol_ylabel);
    set(gca, 'FontSize', fontAxes);

    % Fig 1B: (RIGHT): Std deviation (DET/DIU/DDU) 
    nexttile; hold on; grid on;

    % DET: zero line
    h1s = plot(tt_all, zeros(size(tt_all)), 'k-', 'LineWidth', 3);
    % DIU
    h2s = plot(tt_all, std_M2_all(:,i), 'LineWidth', 3, 'Color', blueBound);
    % DDU
    h3s = plot(tt_all, std_M3_all(:,i), 'LineWidth', 3, 'Color', redBound);

    if i == 1
        hM1_std = h1s;
        hM2_std  = h2s;
        hM3_std  = h3s;
    else
        set(h1s, 'HandleVisibility','off');
        set(h2s, 'HandleVisibility','off');
        set(h3s, 'HandleVisibility','off');
    end

    title(sprintf('Unit %d Inflow Std. Dev.', unitIdx(i)), 'FontSize', fontTitle);
    xlabel('Time (h)');
    ylabel('\sigma_t');
    xlim([1, T])
    ylim([-.001 0.03]);
    set(gca, 'FontSize', fontAxes);
end

sgtitle('State-of-Charge Bounds and Inflow Uncertainty', ...
        'FontSize', fontTitle+10);

% Legends attached to the FIGURE (not the last axes)
% Bounds legend (left side)
lgBounds = legend([hM1_traj hM2_traj hM3_traj], ...
    {'DET', 'DIU', 'DDU'}, ...
    'Orientation','horizontal', 'FontSize', fontLegend, 'Box','on');
drawnow;
posB = lgBounds.Position;
lgBounds.Position = [0.25 - posB(3)/2, 0.005, posB(3), posB(4)];

% Std legend (right side)
lgStd = legend([hM1_std hM2_std hM3_std], ...
    {'DET', 'DIU', 'DDU'}, ...
    'Orientation','horizontal', 'FontSize', fontLegend, 'Box','on');
drawnow;
posS = lgStd.Position;
lgStd.Position = [0.75 - posS(3)/2, 0.005, posS(3), posS(4)];

% saveas(fCombined, fullfile(resultsPath, 'allUnits_bounds_and_std_diu_vs_ddu.png'));

dt = 1;  % hours per time step (change if different)

E_M1 = cumsum(P_M1_allUnits) * dt;
E_M2 = cumsum(P_M2_allUnits) * dt;
E_M3 = cumsum(P_M3_allUnits) * dt;

E_all = [E_M1, E_M2, E_M3];

%{
f4 = figure('Position',[150 300 1100 480]); hold on; grid on;

% Fig 2A: ACCUMULATED ENERGY (Bar Plot)
%{
bE = bar(tt_all, E_all, 'grouped');
bE(1).FaceColor = blackBound;
bE(2).FaceColor = blueBound;
bE(3).FaceColor = redBound;
%}


% Fig 2B: ACCUMULATED ENERGY (Shaded area line Plot)

% Shaded areas
a1 = area(tt_all, E_M1);
a2 = area(tt_all, E_M2);
a3 = area(tt_all, E_M3);

a1.FaceColor = blackBound;
a2.FaceColor = blueBound;
a3.FaceColor = redBound;

a1.FaceAlpha = 0.4;
a2.FaceAlpha = 0.15;
a3.FaceAlpha = 0.15;

a1.EdgeColor = 'none';
a2.EdgeColor = 'none';
a3.EdgeColor = 'none';

l1 = plot(tt_all, E_M1, 'LineWidth', 2, 'Color', blackBound);
l2 = plot(tt_all, E_M2, 'LineWidth', 2, 'Color', blueBound);
l3 = plot(tt_all, E_M3, 'LineWidth', 2, 'Color', redBound);

ylabel('Accumulated Energy (p.u.·h)');
xlabel('Time (h)');
xlim([1 T])
set(gca, 'FontSize', fontAxes);
set(findall(gcf,'Type','text'), 'FontSize', fontAxes);
legend([l1 l2 l3], {'DET', 'DIU', 'DDU'}, ...
       'Location', 'best', 'FontSize', fontLegend);
title('Accumulated Energy Dispatch', 'FontSize', fontTitle);

% saveas(f4, fullfile(resultsPath, 'allUnits_accumulated_energy.png'));
%}
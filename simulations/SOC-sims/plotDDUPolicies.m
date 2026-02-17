%% Plot Driver for SSH vs BON under DDU Comparison
% ========================================================================
% Author: Eliza Cohn
% Description: Overlays policy behavior derived during optimization sims
%              for all units. For each unit k, loads that unit's SSH and BON
%              optimization results (DDU only) and produces:
%                 1) Release
%                 2) Hydraulic Head
% ========================================================================

close all; clc;

% Loading parameters
pathSSH = "./resultsSSH/wet";
pathBON = "./resultsBonferroni/wet";
tag     = 'ddu';      % compare under DDU only

font  = 16;

% Load files for Unit 1 just to get sysparams / n_units / T
D_SSH = load(fullfile(pathSSH, sprintf('results_unit1_%s.mat', lower(tag))));
D_BON = load(fullfile(pathBON, sprintf('results_unit1_%s.mat', lower(tag))));

n_units = numel(D_SSH.sysparams);
T       = D_SSH.T;

%% Figure 1: Comparing SSH vs BON (DDU)
fig = figure('Position',[100 100 1600 300*n_units]);

subfig_n = 2;  % number of states to plot

sgtitle('Optimal Policy Trajectories (DDU): SSH vs BON', ...
    'FontSize', font+6, 'FontWeight','bold');

labels = {'SSH','BON'};

% For legend: store handles from the first unit's u-plot
hSSH_leg = []; hBON_leg = [];

for i = 1:n_units
    sp   = D_SSH.sysparams(i);
    base = (i-1)*5;

    % Extract trajectory from SSH (DDU)
    V_SSH = D_SSH.X(:, base+1);
    u_SSH = D_SSH.X(:, base+3);

    % Extract trajectory from BON (DDU)
    V_BON = D_BON.X(:, base+1);
    u_BON = D_BON.X(:, base+3);

    % Calculate hydraulic head
    headSSH = sp.a .* (V_SSH.^sp.b);
    headBON = sp.a .* (V_BON.^sp.b);

    % Row index for subplots (1..n_units)
    row = i - 1;

    % u_t
    ax = subplot(n_units, subfig_n, row*subfig_n + 1);
    hSSH = plot(u_SSH,'-r','LineWidth',2); hold on;
    hBON = plot(u_BON,'-b','LineWidth',2);
    yline(sp.max_ut,'--k'); yline(sp.min_ut,'--k');
    if i == 1
        title('Water Release [m^3]','FontSize',font);
    end
    if i == n_units
        xlabel('Time [hr]');
    end
    ylabel(sprintf('Unit %d', sp.unit),'FontWeight','bold');
    xlim([1, T]);
    set(ax,'FontSize',font);

    % store legend handles from first unit only
    if i == 1
        hSSH_leg = hSSH; 
        hBON_leg = hBON;
    end

    % head_t
    ax = subplot(n_units, subfig_n, row*subfig_n + 2);
    plot(headSSH,'-r','LineWidth',2); hold on;
    plot(headBON,'-b','LineWidth',2);
    yline(sp.max_h,'--k'); yline(sp.min_h,'--k');
    if i == 1
        title('Hydraulic Head [m]','FontSize',font);
    end
    if i == n_units
        xlabel('Time [hr]');
    end
    xlim([1, T]);
    set(ax,'FontSize',font);
end

% Global legend
if ~isempty(hSSH_leg)
    lg = legend([hSSH_leg hBON_leg], labels{:}, ...
                'Orientation', 'horizontal', ...
                'FontSize', 10, ...
                'Box', 'on');

    % Center the legend dynamically (same logic you used)
    drawnow;
    legendPos = lg.Position;
    legendWidth = legendPos(3);
    lg.Position = [0.5 - legendWidth/2, 0.01, legendWidth, legendPos(4)];
end


%% Summary Statistics 

algNames = {'SSH','BON'};
Data     = {D_SSH, D_BON};

fprintf('DDU POLICY COMPARISON: SSH vs BON\n');

for a = 1:2
    D = Data{a};

    total_gen   = 0;
    total_water = 0;

    fprintf('\n%s:\n', algNames{a});

    for i = 1:n_units
        base = (i-1)*5;

        gen_i = sum(D.X(:, base+2));   % total generation
        wat_i = sum(D.X(:, base+3));   % total release

        total_gen   = total_gen   + gen_i;
        total_water = total_water + wat_i;

        fprintf('Unit %d Generation: %.4f\n', i, gen_i);
    end

    efficiency = total_gen / total_water;

    fprintf('Total Generation: %.4f\n', total_gen);
    fprintf('Total Release: %.4f\n', total_water);
    fprintf('System Efficiency: %.4f\n', efficiency);
end


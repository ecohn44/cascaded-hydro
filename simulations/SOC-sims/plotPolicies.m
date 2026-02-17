%% Plot Driver for DET vs DIU vs DDU Comparison 
% ========================================================================
% Author: Eliza Cohn
% Description: Overlays policy behavior derived during optimization sims
%              for all units. For each unit k, loads that unit's DET, DIU,
%              and DDU optimization results separately and produces:
%                 1) Release
%                 2) Volume
% ========================================================================

close all; clc;

season = "wet";

% Loading parameters 
path = "./resultsSSH/" + season;
tag1 = 'det'; tag2 = 'diu'; tag3 = 'ddu';
printplot = false;

font  = 24;

% Load files for Unit 1 just to get sysparams / n_units
D1 = load(fullfile(path, sprintf('results_unit1_%s.mat', lower(tag1))));
D2 = load(fullfile(path, sprintf('results_unit1_%s.mat', lower(tag2))));
D3 = load(fullfile(path, sprintf('results_unit1_%s.mat', lower(tag3))));

% Define colors 
c1 = [0.0000, 0.4470, 0.7410]; 
c2 = [0.6350, 0.0780, 0.1840];
c3 = [0.4660, 0.6740, 0.1880];

n_units = numel(D1.sysparams);
T       = D1.T;

%% Figure 1: Comparing Uncertainty Frameworks 
% Create one big figure for all units
fig = figure('Position',[100 100 400*n_units 300*n_units]);

if season == "dry"
    subfig_n = 2;  % number of states to plot
else
    subfig_n = 3; 
end

% 'Optimal Policy Trajectories under Uncertainty Frameworks'
if path == "./resultsBonferroni/"+season
    alg = "Bonferroni";
    sgtitle('Bonferroni Trajectories', ...
        'FontSize', 36, 'FontWeight','bold');
else
    alg = "SSH";
    sgtitle('Supporting Hyperplane Trajectories', ...
        'FontSize', 36, 'FontWeight','bold');
end

% Print metrics
printBenchLatexBlock(D1, D2, D3, alg, 5, season);

labels = {upper(tag1), upper(tag2), upper(tag3)};

% For legend: store handles from the first unit's u-plot
h1_leg = []; h2_leg = []; h3_leg = [];

for i = 1:n_units
    sp   = D1.sysparams(i);
    base = (i-1)*5;

    % Extract trajectory from Policy 1
    V1 = D1.X(:, base+1);  p1 = D1.X(:, base+2);
    u1 = D1.X(:, base+3);  s1 = D1.X(:, base+4); 
    V1min = D1.V_eff(:, 2*i);
    V1max = D1.V_eff(:, 2*i - 1);

    % Extract trajectory from Policy 2
    V2 = D2.X(:, base+1);  p2 = D2.X(:, base+2);
    u2 = D2.X(:, base+3);  s2 = D2.X(:, base+4); 
    V2min = D2.V_eff(:, 2*i);
    V2max = D2.V_eff(:, 2*i - 1);

    % Extract trajectory from Policy 3
    V3 = D3.X(:, base+1);  p3 = D3.X(:, base+2);
    u3 = D3.X(:, base+3);  s3 = D3.X(:, base+4);  
    V3min = D3.V_eff(:, 2*i);
    V3max = D3.V_eff(:, 2*i - 1);

    % Calculate hydraulic head
    head1 = sp.a .* (V1.^sp.b);
    head1min = sp.a .* (V1min.^sp.b);
    head1max = sp.a .* (V1max.^sp.b);

    head2 = sp.a .* (V2.^sp.b);
    head2min = sp.a .* (V2min.^sp.b);
    head2max = sp.a .* (V2max.^sp.b);

    head3 = sp.a .* (V3.^sp.b);
    head3min = sp.a .* (V3min.^sp.b);
    head3max = sp.a .* (V3max.^sp.b);

    % Row index for subplots (1..n_units)
    row = i - 1;
    % Column mapping:
    % 1: u_t, 2: p_t, 3: s_t, 4: V_t, 5: head_t

    % u_t
    ax = subplot(n_units, subfig_n, row*subfig_n + 1);
    h1 = plot(u1, 'Color', c1, 'LineWidth',2); hold on;
    h2 = plot(u2, 'Color', c2, 'LineWidth',2);
    h3 = plot(u3, 'Color', c3, 'LineWidth',2);
    yline(sp.max_ut,'--k', 'LineWidth',2); yline(sp.min_ut,'--k', 'LineWidth',2);
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
        h1_leg = h1; h2_leg = h2; h3_leg = h3;
    end


    % h_t 
    ax = subplot(n_units, subfig_n, row*subfig_n + 2);
    % yyaxis left
    ax.YColor = 'k'; 
    plot(head1, 'LineStyle','-', 'Color', c1, 'LineWidth',2); hold on;
    plot(head1min, 'LineStyle','--', 'Color', c1, 'LineWidth',2); 
    plot(head1max, 'LineStyle','--', 'Color', c1, 'LineWidth',2); 
    
    plot(head2, 'LineStyle','-', 'Color', c2, 'LineWidth',2);
    plot(head2min, 'LineStyle','--', 'Color', c2, 'LineWidth',2); 
    plot(head2max, 'LineStyle','--', 'Color', c2, 'LineWidth',2);
    
    plot(head3, 'LineStyle','-', 'Color', c3, 'LineWidth',2);
    plot(head3min, 'LineStyle','--', 'Color', c3, 'LineWidth',2); 
    plot(head3max, 'LineStyle','--', 'Color', c3, 'LineWidth',2);
    
    % yline(sp.max_h,'--k'); yline(sp.min_h,'--k');
    
    if i == 1
        title('Hydraulic Head [m]','FontSize',font);
    end 
    if i == n_units
        xlabel('Time [hr]');
    end
    xlim([1, T]);
    if season == "wet"
        ylim([4.3, 5.1]);
    else
        ylim([-.1, 2.5]);
    end 
    set(ax,'FontSize',font);


    if season == "wet"
        ax = subplot(n_units, subfig_n, row*subfig_n + 3);
        plot(s1, 'Color', c1, 'LineWidth',2); hold on;
        plot(s2, 'Color', c2, 'LineWidth',2);
        plot(s3, 'Color', c3, 'LineWidth',2);
        if i == 1
            title('Spill Release [m^3]','FontSize',font);
        end 
        if i == n_units
            xlabel('Time [hr]');
        end
        ylabel(sprintf('Unit %d', sp.unit),'FontWeight','bold');
        xlim([1, T]);
        set(ax,'FontSize',font);
    end


  
end

% Global legend 
if ~isempty(h1_leg)
    lg = legend([h1_leg h2_leg h3_leg], labels{:}, ...
                'Orientation', 'horizontal', ...
                'FontSize', font-2, ...
                'Box', 'on');

    lg.ItemTokenSize = [40, 28];

    % Center the legend dynamically
    drawnow;  % force legend to compute its size
    legendPos = lg.Position;
    legendWidth = legendPos(3);
    lg.Position = [0.5 - legendWidth/2, 0, legendWidth, legendPos(4)];
end

% Save figure if requested
if printplot
    outFile = fullfile("./plots", ...
        sprintf('AllUnits_%s_%s_%s.png', lower(tag1), lower(tag2), lower(tag3)));
    saveas(fig, outFile);
end


function printBenchLatexBlock(Ddet, Ddiu, Dddu, algName, indentMM, season)
% Prints LaTeX block:
% \multirow{...}{*}{Alg}
%  & System Total [MWh]   & DET & DIU & DDU \\
%  & \hspace{Xmm}Unit 1   & ... \\
%  ...
%  & System Spill [m^3]   & ... \\        (wet only)
%  & \hspace{Xmm}Unit 1   & ... \\        (wet only)
%  ...
%  & System Efficiency [MWh/m$^3$] & ... \\
% \hline

    indent = sprintf('\\hspace{%dmm}', indentMM);

    isWet = (season == "wet");

    % Compute [totalE, eff, unitE, totalS, unitS] for each framework
    [E_det, eff_det, uE_det, S_det, uS_det] = metrics(Ddet, isWet);
    [E_diu, eff_diu, uE_diu, S_diu, uS_diu] = metrics(Ddiu, isWet);
    [E_ddu, eff_ddu, uE_ddu, S_ddu, uS_ddu] = metrics(Dddu, isWet);

    n_units = numel(Ddet.sysparams);

    % Multirow count: energy rows (1 + n_units) + (spill rows if wet) + efficiency row
    n_rows = (1 + n_units) + (isWet*(1 + n_units)) + 1;

    fprintf('\\multirow{%d}{*}{%s}\n', n_rows, algName);

    % Energy block 
    fprintf(' & System Total [MWh] & %.2f & %.2f & %.2f \\\\\n', E_det, E_diu, E_ddu);
    for i = 1:n_units
        fprintf(' & %sUnit %d & %.2f & %.2f & %.2f \\\\\n', ...
            indent, i, uE_det(i), uE_diu(i), uE_ddu(i));
    end

    % Spill block (wet season only)
    if isWet
        fprintf(' & System Spill [m$^3$] & %.2f & %.2f & %.2f \\\\\n', S_det, S_diu, S_ddu);
        for i = 1:n_units
            fprintf(' & %sUnit %d & %.2f & %.2f & %.2f \\\\\n', ...
                indent, i, uS_det(i), uS_diu(i), uS_ddu(i));
        end
    end

    % --- Efficiency block ---
    fprintf(' & System Efficiency [MWh/m$^3$] & %.2f & %.2f & %.2f \\\\\n', ...
        eff_det, eff_diu, eff_ddu);

    fprintf('\\hline\n');

    function [totalE, eff, unitE, totalS, unitS] = metrics(D, wantSpill)
        T = D.T;
        n = numel(D.sysparams);

        unitE = zeros(n,1);
        unitS = zeros(n,1);

        totalW = 0;   % total water release (u)
        totalS = 0;   % total spill (s)

        for k = 1:n
            base = (k-1)*5;
            p = D.X(1:T, base+2);   % MWh per hour-step
            u = D.X(1:T, base+3);   % m^3 per hour-step

            unitE(k) = sum(p);
            totalW   = totalW + sum(u);

            if wantSpill
                s = D.X(1:T, base+4);   % m^3 per hour-step
                unitS(k) = sum(s);
                totalS   = totalS + unitS(k);
            end
        end

        totalE = sum(unitE);
        eff    = totalE / totalW;
    end
end

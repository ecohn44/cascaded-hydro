%% Plot Driver for Efficiency Comparison 
% ========================================================================
% Author: Eliza Cohn
% Description: Overlays effiency derived under DIU and DDU (gamma sweep)
% ========================================================================

clc; clear; close all;

effPath = './resultsBonferroni/eff';

% Figure and font settings
fontAxes   = 14;
fontLegend = 13;
fontTitle  = 16;

n_units = 3;

% Aggregate DIU across units (single curve)
P_M2 = [];
U_M2 = [];

% Aggregate DDU across units (one curve per gamma)
gammaList = [];      % unique gammas discovered
P_DDU = {};          % cell per gamma: P sum
U_DDU = {};          % cell per gamma: U sum

%% Simulation Data
for uu = 1:n_units
    
    fprintf('\n Processing Unit %d \n', uu);

    % Load DIU (single file per unit)
    diuFile = fullfile(effPath, sprintf('results_unit%d_diu.mat', uu));
    M2 = load(diuFile);   % DIU

    % Basic time/index info
    T  = M2.T;
    tt = (1:T)';

    if isempty(P_M2)
        P_M2 = zeros(T,1);
        U_M2 = zeros(T,1);
    end

    % X columns: [V, p, u, s, q] per unit
    baseX = (uu-1)*5;
    colP  = baseX + 2;
    colU  = baseX + 3;

    % Add DIU power + release
    P_M2 = P_M2 + M2.X(:,colP);
    U_M2 = U_M2 + M2.X(:,colU);

    % Load ALL DDU gamma runs for this unit (from /eff)
    dduFiles = dir(fullfile(effPath, sprintf('results_unit%d_ddu_*.mat', uu)));

    for ff = 1:length(dduFiles)

        fname = dduFiles(ff).name;
        fpath = fullfile(dduFiles(ff).folder, fname);

        M3 = load(fpath);   % DDU gamma run

        % --- extract gamma ---
        % Preference: variable inside the .mat (if you saved it)
        if isfield(M3,'gamma')
            gamma = M3.gamma;
        else
            % Fallback: parse from filename like results_unit1_ddu_003.mat
            tok = regexp(fname,'ddu_([0-9]+)\.mat','tokens','once');

            gammaStr = tok{1};                 % e.g. '003'
            gammaInt = str2double(gammaStr);   % -> 3
            nDigits  = length(gammaStr);       % -> 3
            
            gamma = gammaInt * 10^(-nDigits);  % -> 0.003
        end

        % If this is a new gamma, initialize its aggregates
        ii = find(gammaList == gamma, 1);
        if isempty(ii)
            gammaList(end+1,1) = gamma; 
            P_DDU{end+1,1} = zeros(T,1); 
            U_DDU{end+1,1} = zeros(T,1);
            ii = length(gammaList);
        end

        % Add this unit's DDU power + release for this gamma
        P_DDU{ii} = P_DDU{ii} + M3.X(:,colP);
        U_DDU{ii} = U_DDU{ii} + M3.X(:,colU);

    end
end

% Sort gammas (so legend + curves are ordered)
[gammaList, ord] = sort(gammaList);
P_DDU = P_DDU(ord);
U_DDU = U_DDU(ord);

gammaList = gammaList(1:end-1);
P_DDU     = P_DDU(1:end-1);
U_DDU     = U_DDU(1:end-1);

%% Figure 1: Policy Effiency (DIU vs DDU(gamma))

k = 9;   % smoothing window (odd integer looks best)

% Efficiency (system-level)
eta_M2   = P_M2 ./ U_M2;
eta_M2_s = smoothdata(eta_M2,'movmean',k);

% Figure setup
figure('Color','w','Position',[200 200 720 360]);
hold on; box on;

% Colors + weights
col_M2 = [0.05 0.30 0.65];    % DIU (blue)
lw = 1.5;

% Plot DIU
plot(tt, eta_M2_s, 'LineWidth', lw, 'Color', col_M2);

% Plot DDU for each gamma
cols = [ ...
    0.70 0.88 0.72;   % very light green
    0.40 0.75 0.45;   % mid green
    0.15 0.55 0.25];  % deep green
leg  = cell(1, 1 + length(gammaList));
leg{1} = 'DIU';

for ii = 1:length(gammaList)

    eta_g   = P_DDU{ii} ./ U_DDU{ii};
    eta_g_s = smoothdata(eta_g,'movmean',k);

    plot(tt, eta_g_s, 'LineWidth', lw, 'Color', cols(ii,:));

    %leg{1+ii} = sprintf('DDU, \\gamma = %g', gammaList(ii));
end

t_drought = 34;

xline(t_drought, '--', ...
      'LineWidth', 1.6, ...
      'Color', [0.85 0.10 0.10]);

% Drought event annotation
text(t_drought + 1, 1.55, 'Drought event', ...
     'Interpreter','latex', ...
     'FontSize', fontLegend, ...
     'Color', [0.85 0.10 0.10], ...
     'HorizontalAlignment','left');

leg = {'DIU','DDU - \gamma_1','DDU - \gamma_2','DDU - \gamma_3'};

% Axes labels
xlabel('Time (hr)','Interpreter','latex','FontSize',fontAxes);
ylabel('System Efficiency (MWh/m$^3$)','Interpreter','latex','FontSize',fontAxes);

% Legend
lgd = legend(leg, 'Location','northeast', 'FontSize',fontLegend, 'Box','on');

pos = lgd.Position;   % [x y w h] in normalized figure units

annotation('textbox', ...
    [pos(1), pos(2)-0.06, pos(3), 0.06], ...
    'String', '$\gamma_1 < \gamma_2 < \gamma_3$', ...
    'Interpreter','latex', ...
    'FontSize',16, ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center');

xlim([1, T]);

% Axes styling
set(gca,'FontSize',fontAxes, 'LineWidth',1.2, 'TickLabelInterpreter','latex');
grid off;

% Export
exportgraphics(gcf,'figures/efficiency.png','Resolution',300);

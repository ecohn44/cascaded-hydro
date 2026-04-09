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

season = "dry";
alg = "SSH";

N = 40;
n = 3;

% Load inflow data 
[modelparams, sysparams, seasonparams] = dataload(n, N);

% Scale Metrics
% Historical anchors
scale.dry_muQ = 3000;
scale.dry_stdQ = 686;
scale.mu_Q = scale.dry_muQ;
scale.std_Q = scale.dry_stdQ;

% Simulated metrics
scale.mu_q = seasonparams.q0;
scale.std_q = 0.004;

% Hydraulic head scales
scale.H0 = 10;
scale.dH_m = 10; 
scale.alpha_q = scale.std_Q / scale.std_q;   
scale.eta = 0.9;

% Loading parameters 
path = "./resultsSSH/" + season;
tag1 = 'det'; tag2 = 'diu'; tag3 = 'ddu';
printplot = false;

set(groot,'defaultAxesFontName','Helvetica');
set(groot,'defaultTextFontName','Helvetica');
font_title = 14;
font_axis  = 14;
font_tick  = 14;
font_unit  = 14;
font_leg   = 14;
lw    = 2.5;

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

% Create one big figure for all units
fig = figure('Position',[100 100 900 700]);

subfig_n = 2;  % number of cols


labels = {upper(tag1), upper(tag2), upper(tag3)};

% Store legend handles 
h1_leg = []; h2_leg = []; h3_leg = [];

% Precompute column-wise limits for consistent y-axes 
u_all = [];  H_all = [];
for i = 1:n_units
    sp   = D1.sysparams(i);
    base = (i-1)*5;

    % Convert releases (normalized) -> physical 10^3 m^3/s
    u1 = D1.X(:, base+3);  u2 = D2.X(:, base+3);  u3 = D3.X(:, base+3);

    u_all = [u_all; (scale.alpha_q*u1)/1e3; (scale.alpha_q*u2)/1e3; (scale.alpha_q*u3)/1e3; ...
                   (scale.alpha_q*sp.min_ut)/1e3; (scale.alpha_q*sp.max_ut)/1e3];
end

% Y-tick formatting (release)
u_max = max(u_all);
u_step = 2;                        
u_hi = ceil(u_max/u_step)*u_step;
u_ticks = 0:u_step:u_hi;

% Y-tick formatting (head)
H_min = 5; H_max = 12;
H_step = 2;                         
H_lo = floor(H_min/H_step)*H_step;
H_hi = ceil(H_max/H_step)*H_step;
H_ticks = H_lo:H_step:H_hi;

for i = 1:n_units
    sp   = D1.sysparams(i);
    base = (i-1)*5;

    % Extract trajectories
    V1 = D1.X(:, base+1);  u1 = D1.X(:, base+3);  s1 = D1.X(:, base+4); 
    V2 = D2.X(:, base+1);  u2 = D2.X(:, base+3);  s2 = D2.X(:, base+4); 
    V3 = D3.X(:, base+1);  u3 = D3.X(:, base+3);  s3 = D3.X(:, base+4);  
    V1min = D1.V_eff(:, 2*i); V2min = D2.V_eff(:, 2*i); V3min = D3.V_eff(:, 2*i);

    % Calculate hydraulic head (model proxy)
    head1 = sp.a .* (V1.^sp.b);
    head1min = sp.a .* (V1min.^sp.b);
    head2 = sp.a .* (V2.^sp.b);
    head2min = sp.a .* (V2min.^sp.b);
    head3 = sp.a .* (V3.^sp.b);
    head3min = sp.a .* (V3min.^sp.b);

    row = i - 1;

    % Release: physical scaling
    u1_phys = (scale.alpha_q * u1) / 1e3; 
    u2_phys = (scale.alpha_q * u2) / 1e3; 
    u3_phys = (scale.alpha_q * u3) / 1e3; 

    ax = subplot(n_units, subfig_n, row*subfig_n + 1);
    h1 = plot(u1_phys, 'Color', c1, 'LineWidth', lw); hold on;
    h2 = plot(u2_phys, 'Color', c2, 'LineWidth', lw);
    h3 = plot(u3_phys, 'Color', c3, 'LineWidth', lw);

    yline((scale.alpha_q*sp.max_ut)/1e3,'--','Color',[0.3 0.3 0.3],'LineWidth',1.5)
    yline((scale.alpha_q*sp.min_ut)/1e3,'--','Color',[0.3 0.3 0.3],'LineWidth',1.5)

    % column header only (small + normal weight)
    if i == 1
        title('Release (10^3 m^3/s)','FontSize',font_title,'FontWeight','normal');
    elseif i == n_units
        xlabel('Time (hours)','FontSize',font_axis);
    end
    ylabel(sprintf('Unit %d', sp.unit),'FontSize',font_unit);

    % xlim([1, T]);
    % ylim([0, u_hi]);
    ylim([5 10])
    xlim([30 50])
    
    set(gca,'FontSize',12,'LineWidth',1,'Box','on','TickDir','out')

    yticks(u_ticks);

    grid on; set(ax,'XGrid','off','YGrid','on');
    set(ax,'FontSize',font_tick);

    if i == 1
        h1_leg = h1; h2_leg = h2; h3_leg = h3;
    end

    % Head: percent-full mapping -> physical meters 
    hmin_model = 0;           
    hmax_model = sp.a;       

    head1_norm = min(max((head1 - hmin_model) ./ (hmax_model - hmin_model), 0), 1);
    head1_minnorm = min(max((head1min - hmin_model) ./ (hmax_model - hmin_model), 0), 1);
    
    head2_norm = min(max((head2 - hmin_model) ./ (hmax_model - hmin_model), 0), 1);
    head2_minnorm = min(max((head2min - hmin_model) ./ (hmax_model - hmin_model), 0), 1);
    
    head3_norm = min(max((head3 - hmin_model) ./ (hmax_model - hmin_model), 0), 1);
    head3_minnorm = min(max((head3min - hmin_model) ./ (hmax_model - hmin_model), 0), 1);

    head1_phys = scale.H0 + (head1_norm - 0.5) * scale.dH_m; 
    head1_minphys = scale.H0 + (head1_minnorm - 0.5) * scale.dH_m; 

    head2_phys = scale.H0 + (head2_norm - 0.5) * scale.dH_m;  
    head2_minphys = scale.H0 + (head2_minnorm - 0.5) * scale.dH_m; 

    head3_phys = scale.H0 + (head3_norm - 0.5) * scale.dH_m; 
    head3_minphys = scale.H0 + (head3_minnorm - 0.5) * scale.dH_m; 

    ax = subplot(n_units, subfig_n, row*subfig_n + 2);
    plot(head1_phys, 'LineStyle','-', 'Color', c1, 'LineWidth', lw); hold on;
    plot(head1_minphys, 'LineStyle','--', 'Color', c1, 'LineWidth', 1.5);
    plot(head2_phys, 'LineStyle','-', 'Color', c2, 'LineWidth', lw);
    plot(head2_minphys, 'LineStyle','--', 'Color', c2, 'LineWidth', 1.5);
    plot(head3_phys, 'LineStyle','-', 'Color', c3, 'LineWidth', lw);
    plot(head3_minphys, 'LineStyle','--', 'Color', c3, 'LineWidth', 1.5);

    if i == 1
        title('Head (m)','FontSize',font_title,'FontWeight','normal');
    end 
    if i == n_units
        xlabel('Time (hours)');
    end

    xlim([1, T]);
    ylim([H_lo, H_hi]);
    yticks(H_ticks);

    grid on; set(ax,'XGrid','off','YGrid','on');
    set(ax,'FontSize',font_tick);

end

% Global legend (below)
if ~isempty(h1_leg)
    lg = legend([h1_leg h2_leg h3_leg], labels{:}, ...
                'Orientation', 'horizontal', ...
                'FontSize', font_leg, ...
                'Box', 'off');
    lg.ItemTokenSize = [40, 28];

    drawnow;
    legendPos = lg.Position;
    legendWidth = legendPos(3);
    lg.Position = [0.5 - legendWidth/2, 0.01, legendWidth, legendPos(4)];
end

set(gcf,'Renderer','painters');
exportgraphics(gcf,'figures/dry_policy_trajectories.pdf','ContentType','vector');


% Global legend (bottom)
labels = {upper(tag1), upper(tag2), upper(tag3)};
lg = legend([h1_leg h2_leg h3_leg], labels{:}, ...
            'Orientation', 'horizontal', ...
            'FontSize', font_leg, ...
            'Box', 'off');
lg.ItemTokenSize = [40, 28];
drawnow;
legendPos = lg.Position;
legendWidth = legendPos(3);
lg.Position = [0.5 - legendWidth/2, 0.01, legendWidth, legendPos(4)];

printBenchLatexBlock(D1, D2, D3, alg, 5, scale)



function printBenchLatexBlock(Ddet, Ddiu, Dddu, algName, indentMM, scale)

    dt_hr = 1;       
    rho   = 1000;    
    g     = 9.81;    
    eta   = scale.eta;

    indent = sprintf('\\hspace{%dmm}', indentMM);

    % Compute energy + released volume for each framework
    [E_det, uE_det, Vrel_det] = energySum(Ddet);
    [E_diu, uE_diu, Vrel_diu] = energySum(Ddiu);
    [E_ddu, uE_ddu, Vrel_ddu] = energySum(Dddu);

    n_units = numel(Ddet.sysparams);
    T_hours = Ddet.T * dt_hr;

    % Average powers
    Pavg_sys_det = E_det / T_hours;
    Pavg_sys_diu = E_diu / T_hours;
    Pavg_sys_ddu = E_ddu / T_hours;

    Pavg_unit_det = Pavg_sys_det / n_units;
    Pavg_unit_diu = Pavg_sys_diu / n_units;
    Pavg_unit_ddu = Pavg_sys_ddu / n_units;

    % System efficiency [MWh/m^3]
    Eff_det = E_det / (Vrel_det / 1e6);
    Eff_diu = E_diu / (Vrel_diu / 1e6);
    Eff_ddu = E_ddu / (Vrel_ddu / 1e6);

    fmtInt = @(x) regexprep(sprintf('%.0f', round(x)), '\B(?=(\d{3})+(?!\d))', ',');

    % Multirow count: energy rows (1+n_units) + 2 avg-power rows + 1 efficiency row
    n_rows = (1 + n_units) + 3;
    
    fprintf('\\multirow{%d}{*}{%s}\n', n_rows, algName);
    
    fprintf(' & System Total [MWh] & %s & %s & %s \\\\\n', ...
        fmtInt(E_det), fmtInt(E_diu), fmtInt(E_ddu));
    
    for i = 1:n_units
        fprintf(' & %sUnit %d & %s & %s & %s \\\\\n', ...
            indent, i, fmtInt(uE_det(i)), fmtInt(uE_diu(i)), fmtInt(uE_ddu(i)));
    end
    
    fprintf(' & Avg System Power [MW] & %s & %s & %s \\\\\n', ...
        fmtInt(Pavg_sys_det), fmtInt(Pavg_sys_diu), fmtInt(Pavg_sys_ddu));
    
    fprintf(' & Avg Unit Power [MW] & %s & %s & %s \\\\\n', ...
        fmtInt(Pavg_unit_det), fmtInt(Pavg_unit_diu), fmtInt(Pavg_unit_ddu));

    fprintf(' & System Efficiency [MWh/$10^6$ m$^3$] & %.2f & %.2f & %.2f \\\\\n', ...
        Eff_det, Eff_diu, Eff_ddu);
    
    fprintf('\\hline\n');

    function [E_system, E_unit, Vrel_system] = energySum(D)
        T = D.T;
        n = numel(D.sysparams);

        E_unit = zeros(n,1);
        Vrel_unit = zeros(n,1);

        for k = 1:n
            sp   = D.sysparams(k);
            base = (k-1)*5;

            % Turbine release (normalized) -> m^3/s
            u_norm = D.X(1:T, base+3);
            u_m3ps = scale.alpha_q * u_norm;

            % Total released volume [m^3]
            Vrel_unit(k) = sum(u_m3ps) * 3600 * dt_hr;

            % Head proxy from normalized storage -> physical meters
            V_norm  = D.X(1:T, base+1);
            h_model = sp.a .* (V_norm.^sp.b);
            h_norm  = min(max(h_model ./ sp.a, 0), 1);

            H_m = scale.H0 + (h_norm - 0.5) * scale.dH_m;

            % Power -> Energy
            P_MW = (rho * g * eta .* u_m3ps .* H_m) / 1e6;
            E_unit(k) = sum(P_MW) * dt_hr;
        end

        E_system = sum(E_unit);
        Vrel_system = sum(Vrel_unit);
    end
end
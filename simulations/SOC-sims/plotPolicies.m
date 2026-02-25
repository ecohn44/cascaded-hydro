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
alg = "BON";

% Scale Metrics
scale.mu_Q = 4600;
scale.std_Q = 2100;
scale.mu_q   = 0.0406; 
scale.std_q  = 0.0060; 
scale.H0 = 10;
scale.dH_m = 10; 
scale.alpha_q = scale.std_Q / scale.std_q;   
scale.eta = 0.9;

% Loading parameters 
path = "./resultsBonferroni/" + season;
tag1 = 'det'; tag2 = 'diu'; tag3 = 'ddu';
printplot = false;

set(groot,'defaultAxesFontName','Helvetica');
set(groot,'defaultTextFontName','Helvetica');
font_title = 22;
font_axis  = 18;
font_tick  = 16;
font_unit  = 18;
font_leg   = 18;
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
fig = figure('Position',[100 100 400*n_units 300*n_units]);

if season == "dry"
    subfig_n = 2;  % number of states to plot
else
    subfig_n = 3; 
end

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

    % Convert heads (normalized proxy) -> percent-full -> physical meters
    V1 = D1.X(:, base+1);  V2 = D2.X(:, base+1);  V3 = D3.X(:, base+1);
    head1 = sp.a .* (V1.^sp.b);
    head2 = sp.a .* (V2.^sp.b);
    head3 = sp.a .* (V3.^sp.b);

    hmin_model = 0;           
    hmax_model = sp.a;

    h1n = min(max((head1 - hmin_model) ./ (hmax_model - hmin_model), 0), 1);
    h2n = min(max((head2 - hmin_model) ./ (hmax_model - hmin_model), 0), 1);
    h3n = min(max((head3 - hmin_model) ./ (hmax_model - hmin_model), 0), 1);

    H1 = scale.H0 + (h1n - 0.5) * scale.dH_m;
    H2 = scale.H0 + (h2n - 0.5) * scale.dH_m;
    H3 = scale.H0 + (h3n - 0.5) * scale.dH_m;

    H_all = [H_all; H1; H2; H3];
end

% Y-tick formatting (release)
u_max = max(u_all);
u_step = 5;                        
u_hi = ceil(u_max/u_step)*u_step;
u_ticks = 0:u_step:u_hi;

% Y-tick formatting (head)
H_min = min(H_all); H_max = max(H_all);
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

    % Calculate hydraulic head (model proxy)
    head1 = sp.a .* (V1.^sp.b);
    head2 = sp.a .* (V2.^sp.b);
    head3 = sp.a .* (V3.^sp.b);

    row = i - 1;

    % Release: physical scaling
    u1_phys = (scale.alpha_q * u1) / 1e3; 
    u2_phys = (scale.alpha_q * u2) / 1e3; 
    u3_phys = (scale.alpha_q * u3) / 1e3; 

    ax = subplot(n_units, subfig_n, row*subfig_n + 1);
    h1 = plot(u1_phys, 'Color', c1, 'LineWidth', lw); hold on;
    h2 = plot(u2_phys, 'Color', c2, 'LineWidth', lw);
    h3 = plot(u3_phys, 'Color', c3, 'LineWidth', lw);

    yline((scale.alpha_q*sp.max_ut)/1e3,'--','Color',[0.3 0.3 0.3],'LineWidth',1.2)
    yline((scale.alpha_q*sp.min_ut)/1e3,'--','Color',[0.3 0.3 0.3],'LineWidth',1.2)

    % column header only (small + normal weight)
    if i == 1
        title('Release (10^3 m^3/s)','FontSize',font_title,'FontWeight','normal');
    end
    if i == n_units
        xlabel('Time (hours)','FontSize',font_axis);
    end
    ylabel(sprintf('Unit %d', sp.unit),'FontSize',font_unit);

    xlim([1, T]);
    ylim([0, u_hi]);
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
    head2_norm = min(max((head2 - hmin_model) ./ (hmax_model - hmin_model), 0), 1);
    head3_norm = min(max((head3 - hmin_model) ./ (hmax_model - hmin_model), 0), 1);

    head1_phys = scale.H0 + (head1_norm - 0.5) * scale.dH_m;  
    head2_phys = scale.H0 + (head2_norm - 0.5) * scale.dH_m;  
    head3_phys = scale.H0 + (head3_norm - 0.5) * scale.dH_m;  

    ax = subplot(n_units, subfig_n, row*subfig_n + 2);
    plot(head1_phys, 'LineStyle','-', 'Color', c1, 'LineWidth', lw); hold on;
    plot(head2_phys, 'LineStyle','-', 'Color', c2, 'LineWidth', lw); 
    plot(head3_phys, 'LineStyle','-', 'Color', c3, 'LineWidth', lw);

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

    if season == "wet"
        ax = subplot(n_units, subfig_n, row*subfig_n + 3);
        plot(s1, 'Color', c1, 'LineWidth', lw); hold on;
        plot(s2, 'Color', c2, 'LineWidth', lw);
        plot(s3, 'Color', c3, 'LineWidth', lw);
        if i == 1
            title('Spill (m^3/s)','FontSize',font,'FontWeight','normal');
        end 
        if i == n_units
            xlabel('Time (hours)','FontSize',font_axis);
        end
        ylabel(sprintf('Unit %d', sp.unit),'FontSize',font_unit);
        xlim([1, T]);
        grid on; set(ax,'XGrid','off','YGrid','on');
        set(ax,'FontSize',font_tick);
    end
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


%% Fig 2: Accumulated Energy 
rho = 1000; g = 9.81; eta = scale.eta;
tt_all = (1:T)';

sysPowerMW = @(D) local_system_power_MW(D, scale, rho, g, eta);

Psys1 = sysPowerMW(D1);   % DET
Psys2 = sysPowerMW(D2);   % DIU
Psys3 = sysPowerMW(D3);   % DDU

Eacc1 = cumsum(Psys1);    % MWh (dt=1 hr)
Eacc2 = cumsum(Psys2);
Eacc3 = cumsum(Psys3);

% Plot
fE = figure('Position',[100 100 300*n_units 300*n_units]); hold on; grid on;

% Shaded areas (accumulated energy)
a1 = area(tt_all, Eacc1); a1.FaceColor = c1; a1.FaceAlpha = 0.22; a1.EdgeColor = 'none';
a2 = area(tt_all, Eacc2); a2.FaceColor = c2; a2.FaceAlpha = 0.14; a2.EdgeColor = 'none';
a3 = area(tt_all, Eacc3); a3.FaceColor = c3; a3.FaceAlpha = 0.14; a3.EdgeColor = 'none';

% Lines on top
l1 = plot(tt_all, Eacc1, 'LineWidth', lw, 'Color', c1);
l2 = plot(tt_all, Eacc2, 'LineWidth', lw, 'Color', c2);
l3 = plot(tt_all, Eacc3, 'LineWidth', lw, 'Color', c3);

xlabel('Time (hours)','FontSize',font_axis);
ylabel('Accumulated Energy (MWh)','FontSize',font_axis);
set(gca,'FontSize',font_tick);
xlim([1 T]);

% Keep legend consistent with your main figure labels
legend([l1 l2 l3], labels{:}, 'Location','northwest', 'FontSize',font_leg, 'Box','off');

set(gcf,'Renderer','painters');
exportgraphics(gcf,'figures/dry_accumulated_energy.pdf','ContentType','vector');

function P_sys_MW = local_system_power_MW(D, scale, rho, g, eta)
    % Returns system total power time series (MW): sum over units, each hour.
    T = D.T;
    n = numel(D.sysparams);

    P_sys_MW = zeros(T,1);

    for k = 1:n
        sp   = D.sysparams(k);
        base = (k-1)*5;

        % Turbine release (normalized) -> m^3/s
        u_norm = D.X(1:T, base+3);
        Q_m3ps = scale.alpha_q * u_norm;

        % Head proxy -> percent-full -> physical meters
        V_norm  = D.X(1:T, base+1);
        h_model = sp.a .* (V_norm.^sp.b);           % in [0, a]
        h_norm  = min(max(h_model ./ sp.a, 0), 1);  % clip to [0,1]
        H_m     = scale.H0 + (h_norm - 0.5) * scale.dH_m;

        % Unit power (MW)
        P_unit_MW = (rho * g * eta .* Q_m3ps .* H_m) / 1e6;

        % Add to system
        P_sys_MW = P_sys_MW + P_unit_MW;
    end
end


printBenchLatexBlock(D1, D2, D3, alg, 5, scale)



function printBenchLatexBlock(Ddet, Ddiu, Dddu, algName, indentMM, scale)
% =========================================================================
% printBenchLatexBlock (SIMPLIFIED: ENERGY + AVG POWER)
%
% Prints LaTeX rows for:
%   - System Total Energy [MWh]  (sum across units)
%   - Unit i Energy [MWh]
%   - Avg System Power [MW]      (= System Total / T_hours)
%   - Avg Unit Power [MW]        (= Avg System / n_units)
%
% Physical energy computed from post-processed turbine release and head:
%   u_norm  -> u_m3ps = scale.alpha_q * u_norm
%   V_norm  -> h_model = a * V^b
%           -> h_norm  = clip(h_model / a, 0, 1)
%          -> H_m     = scale.H0 + (h_norm - 0.5) * scale.dH_m
%   P_MW(t) = rho*g*eta*u_m3ps(t)*H_m(t) / 1e6
%   E_MWh   = sum_t P_MW(t)   (dt = 1 hour)
% =========================================================================

    dt_hr = 1;       % fixed 1-hour timestep
    rho   = 1000;    % kg/m^3
    g     = 9.81;    % m/s^2
    eta   = scale.eta;

    indent = sprintf('\\hspace{%dmm}', indentMM);

    % Compute energy for each framework
    [E_det, uE_det] = energySum(Ddet);
    [E_diu, uE_diu] = energySum(Ddiu);
    [E_ddu, uE_ddu] = energySum(Dddu);

    n_units = numel(Ddet.sysparams);
    T_hours = Ddet.T * dt_hr;

    % Average powers
    Pavg_sys_det = E_det / T_hours;
    Pavg_sys_diu = E_diu / T_hours;
    Pavg_sys_ddu = E_ddu / T_hours;

    Pavg_unit_det = Pavg_sys_det / n_units;
    Pavg_unit_diu = Pavg_sys_diu / n_units;
    Pavg_unit_ddu = Pavg_sys_ddu / n_units;

    fmtInt = @(x) regexprep(sprintf('%.0f', round(x)), '\B(?=(\d{3})+(?!\d))', ',');

    % Multirow count: energy rows (1+n_units) + 2 avg-power rows
    n_rows = (1 + n_units) + 2;
    
    fprintf('\\multirow{%d}{*}{%s}\n', n_rows, algName);
    
    % System total energy (integer MWh)
    fprintf(' & System Total [MWh] & %s & %s & %s \\\\\n', ...
        fmtInt(E_det), fmtInt(E_diu), fmtInt(E_ddu));
    
    % Unit-level energy (integer MWh)
    for i = 1:n_units
        fprintf(' & %sUnit %d & %s & %s & %s \\\\\n', ...
            indent, i, fmtInt(uE_det(i)), fmtInt(uE_diu(i)), fmtInt(uE_ddu(i)));
    end
    
    % Avg powers (integer MW)
    fprintf(' & Avg System Power [MW] & %s & %s & %s \\\\\n', ...
        fmtInt(Pavg_sys_det), fmtInt(Pavg_sys_diu), fmtInt(Pavg_sys_ddu));
    
    fprintf(' & Avg Unit Power [MW] & %s & %s & %s \\\\\n', ...
        fmtInt(Pavg_unit_det), fmtInt(Pavg_unit_diu), fmtInt(Pavg_unit_ddu));
    
    fprintf('\\hline\n');

    function [E_system, E_unit] = energySum(D)
        T = D.T;
        n = numel(D.sysparams);

        E_unit = zeros(n,1);

        for k = 1:n
            sp   = D.sysparams(k);
            base = (k-1)*5;

            % Turbine release (normalized) -> m^3/s
            u_norm = D.X(1:T, base+3);
            u_m3ps = scale.alpha_q * u_norm;

            % Head proxy from normalized storage -> percent-full -> physical meters
            V_norm  = D.X(1:T, base+1);
            h_model = sp.a .* (V_norm.^sp.b);                % in [0, a]
            h_norm  = min(max(h_model ./ sp.a, 0), 1);       % clip to [0,1]

            H_m = scale.H0 + (h_norm - 0.5) * scale.dH_m;

            % Instantaneous power in MW, then energy in MWh (dt_hr = 1)
            P_MW = (rho * g * eta .* u_m3ps .* H_m) / 1e6;
            E_unit(k) = sum(P_MW) * dt_hr;   % MWh
        end

        E_system = sum(E_unit);
    end
end
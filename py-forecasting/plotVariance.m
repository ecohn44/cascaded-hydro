%% Author: Eliza Cohn
% Date: February 2026
% Description: Compute residual variance vs upstream discharge proxy
% Paper: EXPLICIT RISK ALLOCATION FOR CASCADED HYDROELECTRIC SYSTEMS UNDER EXTREME EVENTS


close all; clear; clc;

bon_path = "/Users/elizacohn/Desktop/cascaded-hydro/streamflow-data-raw/bonneville/";
tda_path = "/Users/elizacohn/Desktop/cascaded-hydro/streamflow-data-raw/dalles/";

bon_file = bon_path + "bon-fullflow.csv";
tda_file = tda_path + "tda-fullflow.csv";

%% STEP 1: Data Load and Preprocessing

bon = readtable(bon_file);
tda = readtable(tda_file);

% Fix timestamps (same as your Python code)
bon.DateTime = datetime(bon.("DateTime"));
tda.DateTime = datetime(tda.("DateTime"));

% Extract exact columns you used in Python
bon_inflow = bon.("BON_Flow_In_Inst__6Hours_0_RFC_FCST_kcfs_");
tda_outflow = tda.("TDA_Flow_Out_Ave_1Hour_1Hour_CBT_REV_kcfs_");

% Convert kcfs → m³/s 
kcfs_to_m3s = @(x) x * 28.32;
bon_inflow = kcfs_to_m3s(bon_inflow);
tda_outflow = kcfs_to_m3s(tda_outflow);

% Merge datasets
flow = innerjoin( ...
    table(bon.DateTime, bon_inflow, 'VariableNames', {'DateTime','inflow_m3s'}), ...
    table(tda.DateTime, tda_outflow, 'VariableNames', {'DateTime','outflow_m3s'}), ...
    'Keys','DateTime');

% Drop NaNs and convert to 6-hour timestep
flow = flow(~isnan(flow.inflow_m3s), :);
flow = flow(~isnan(flow.outflow_m3s), :);

% Extract aligned signals
t      = flow.DateTime;
Q_down = flow.inflow_m3s;
Q_up   = flow.outflow_m3s;
tau    = 1;        % Use 6-hour lag = 1 step

%% Step 2: Build AR Model 

% Define lagged upstream release predictor (Qup_{t-1})
Q_proxy = [NaN; Q_up(1:end-1)];

% Define lagged AR predictor for valid data (Qdn_{t-1})
Q_down_lag = [NaN; Q_down(1:end-1)];
valid = ~isnan(Q_down_lag);

% Perform the OLS regression and solve for coefficients 
beta = [ones(sum(valid),1), Q_down_lag(valid)] \ Q_down(valid);

% Calculated predicted inflow from fitted OLS model 
Q_pred = beta(1) + beta(2)*Q_down_lag;

% Calculate residuals 
res  = Q_down - Q_pred;
res2 = res.^2;

% Remove final NaN from lag structures and align signals 
valid = ~isnan(Q_proxy) & ~isnan(res2);
Q_proxy = Q_proxy(valid);
res2    = res2(valid);


%% Plot #1: Scattered Residuals vs Upstream Release Magnitude 
figure('Color','w'); hold on

scatter(Q_proxy, res(valid), 8, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);

xlabel('Upstream Release (m^3/s)')
ylabel('Variance ((m^3/s)^2)')
title('Residual Variance vs Upstream Release')

grid on
box on
set(gca,'FontSize',14)


%% Plot #2: Binned Variance vs Upstream Release Magnitude 

nbins = 15;
edges = quantile(Q_proxy, linspace(0,1,nbins+1));

bin_centers = zeros(nbins,1);
var_binned  = zeros(nbins,1);

for i=1:nbins
    
    idx = Q_proxy>=edges(i) & Q_proxy<edges(i+1);
    
    bin_centers(i) = mean(Q_proxy(idx));
    var_binned(i)  = mean(res2(idx));
    
end

figure('Color','w'); hold on

plot(bin_centers, var_binned, 'o-', ...
    'LineWidth',2,...
    'MarkerSize',8)

xlabel('Upstream Release (m^3/s)')
ylabel('Variance ((m^3/s)^2)')
title('Binned Residual Variance vs Upstream Release')

grid on
box on
set(gca,'FontSize',14)

%% Fig 3: Quantile-binned box-bars for Residuals Standard Deviation vs Upstream Release

res_std = abs(res); 

% Quantile bins (equal counts) 
nbins = 12;  
edges = quantile(Q_proxy, linspace(0,1,nbins+1));
edges(1)   = -inf;     % include min safely
edges(end) =  inf;     % include max safely

bin_id = discretize(Q_proxy, edges);
ok = ~isnan(bin_id) & ~isnan(res_std);
bin_id = bin_id(ok);
res_ok = res_std(ok);
Q_ok   = Q_proxy(ok);

% Evenly spaced x positions
xpos = 1:nbins;

figure('Color','w'); clf; hold on

% Visual widths in "index space"
w   = 0.33;    % half-width of box
cap = 0.12;    % half-width of whisker caps

for i = 1:nbins
    r = res_ok(bin_id == i);
    if isempty(r), continue; end

    % Winsorize
    loW = prctile(r, 5);
    hiW = prctile(r, 95);
    r   = min(max(r, loW), hiW);

    % Stats: box = IQR, whiskers = 10–90, median line
    p15 = prctile(r,15);
    p25 = prctile(r,25);
    p50 = prctile(r,50);
    p75 = prctile(r,75);
    p85 = prctile(r,85);

    x = xpos(i);

    fill([x-w x+w x+w x-w], [p25 p25 p75 p75], ...
         [0.35 0.75 1.0], 'FaceAlpha',0.60, ...
         'EdgeColor','k','LineWidth',1.2);

    plot([x-w x+w], [p50 p50], 'r', 'LineWidth',2.2);

    plot([x x], [p15 p25], 'k', 'LineWidth',1.2);
    plot([x x], [p75 p85], 'k', 'LineWidth',1.2);

    plot([x-cap x+cap], [p15 p15], 'k', 'LineWidth',1.2);
    plot([x-cap x+cap], [p85 p85], 'k', 'LineWidth',1.2);
end

grid on; box on
set(gca,'FontSize',14,'LineWidth',1.2)

% Create Legend
h_box = fill(nan, nan, [0.35 0.75 1.0], ...
    'FaceAlpha',0.60, ...
    'EdgeColor','k', ...
    'LineWidth',1.2, ...
    'DisplayName','IQR (25–75%)');

h_med = plot(nan, nan, 'r', ...
    'LineWidth',2.2, ...
    'DisplayName','Median');

h_whisk = plot(nan, nan, 'k', ...
    'LineWidth',1.2, ...
    'DisplayName','Whiskers (15–85%)');

legend([h_box, h_med, h_whisk], 'Location','northwest');

xlabel('Upstream Release (10^3 m^3/s)')  % median bin
ylabel('Residual Magnitude (m^3/s)')
% title('Residual Standard Deviation vs Binned Upstream Release') % quantile bins

xlim([0.5 nbins+0.5])
xticks(xpos)

cent = zeros(1,nbins);
for i=1:nbins
    cent(i) = median(Q_ok(bin_id==i));
end
xticklabels(compose('%.1f', cent/1e3)) 


%% Plot #4: Normal Distributions from DIU and DDU

e = res;   % signed residuals (m^3/s)

% Quantile bins for low/high release regimes
nbins = 12;
qedges = quantile(Q_proxy, linspace(0,1,nbins+1));
qedges(1)   = -inf;
qedges(end) =  inf;

bin_id = discretize(Q_proxy, qedges);

ok = ~isnan(bin_id) & ~isnan(e) & ~isnan(Q_proxy);
bin_id = bin_id(ok);
e_ok   = e(ok);

lowBin  = 1;
highBin = nbins;

% Sigma estimates
sig_low  = std(e_ok(bin_id == lowBin),  0, 'omitnan');   % DDU @ low release bin
sig_high = std(e_ok(bin_id == highBin), 0, 'omitnan');   % DDU @ high release bin
sig_diu  = std(e_ok,                    0, 'omitnan');   % DIU (static)

% Guard against degenerate cases
sig_low  = max(sig_low,  eps);
sig_high = max(sig_high, eps);
sig_diu  = max(sig_diu,  eps);

% X-grid for PDFs
sig_max = max([sig_low, sig_high, sig_diu]);
x = linspace(-4*sig_max, 4*sig_max, 600);

% Normal PDFs (mean 0)
pdf_low  = (1/(sig_low  *sqrt(2*pi))) * exp(-0.5*(x./sig_low ).^2);
pdf_high = (1/(sig_high *sqrt(2*pi))) * exp(-0.5*(x./sig_high).^2);
pdf_diu  = (1/(sig_diu  *sqrt(2*pi))) * exp(-0.5*(x./sig_diu ).^2);

% Color
blue   = [0.0000, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];

% Plot
figure('Color','w'); hold on
plot(x, pdf_low,  'Color', orange,  'LineWidth',2.5);                 % DDU low release
plot(x, pdf_high, 'Color', blue,  'LineWidth',2.5);     % DDU high release
plot(x, pdf_diu,  'k', 'LineWidth',2.5);            % DIU (static)

xlabel('Forecast Error (m$^3$/s)','Interpreter','latex')
ylabel('PDF','Interpreter','latex')
grid on; box on
set(gca,'FontSize',14,'LineWidth',1.2)

% Bin medians for legend context (optional but nice)
Q_low_med  = median(Q_proxy(ok & bin_id==lowBin),  'omitnan');   % careful: bin_id already filtered
Q_high_med = median(Q_proxy(ok & bin_id==highBin), 'omitnan');

legend({ ...
    sprintf('DDU (Low Release)'), ...
    sprintf('DDU (High Release)'), ...
    sprintf('DIU')}, ...
    'Location','northeast');

title('');
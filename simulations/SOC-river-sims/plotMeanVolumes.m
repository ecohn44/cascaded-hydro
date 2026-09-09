%% Plot Driver for DIU vs DDU: Average Across Years
clear; clc; close all;

% User settings
results_file = fullfile("resultsBonferroni","monteCarloResultskrangeMC.mat");
kappa_value = 1;
theta_value = 1;

font = 16;

% Load results
S = load(results_file,"results");
results = S.results;

k = find(results.kappa == kappa_value,1);
h = find(results.thetas == theta_value,1);
m_diu = find(results.frameworks == "diu",1);
m_ddu = find(results.frameworks == "ddu",1);

T = size(results.V,2);
n_units = size(results.V,1);
unitNames = ["McNary","John Day","The Dalles","Bonneville"];
x = 1:T;

% Colors
diu_color  = [0.00 0.45 0.74];
ddu_color  = [0.85 0.33 0.10];
band_color = [0.80 0.80 0.80];

% Plot
fig = figure("Color","w","Position",[100 100 1200 900]);
tl = tiledlayout(n_units,2,"TileSpacing","compact","Padding","compact");

for i = 1:n_units
    sp = results.sysparams(i);

    V_mean = squeeze(mean(results.SOC_mean(i,:,:),3));
    V_p10  = squeeze(mean(results.SOC_p10(i,:,:),3));
    V_p90  = squeeze(mean(results.SOC_p90(i,:,:),3));

    V_diu = squeeze(mean(results.V(i,:,: ,h,m_diu,k,:),[3 7]));
    V_ddu = squeeze(mean(results.V(i,:,: ,h,m_ddu,k,:),[3 7]));

    u_diu = squeeze(mean(results.u(i,:,: ,h,m_diu,k,:),[3 7]));
    u_ddu = squeeze(mean(results.u(i,:,: ,h,m_ddu,k,:),[3 7]));

    head_mean = sp.a .* V_mean.^sp.b;
    head_p10  = sp.a .* V_p10.^sp.b;
    head_p90  = sp.a .* V_p90.^sp.b;
    head_diu  = sp.a .* V_diu.^sp.b;
    head_ddu  = sp.a .* V_ddu.^sp.b;

    % Release
    ax = nexttile((i-1)*2 + 1);
    p_diu = plot(x,u_diu,"Color",diu_color,"LineWidth",3.5); hold on;
    p_ddu = plot(x,u_ddu,"Color",ddu_color,"LineWidth",3.5);
    yline(sp.max_ut,"--k","LineWidth",1.2);
    yline(sp.min_ut,"--k","LineWidth",1.2);

    if i == 1
        title("Release (m^3/s)");
        h_diu_leg = p_diu;
        h_ddu_leg = p_ddu;
    end
    if i == n_units
        xlabel("Time (hr)");
    end

    ylabel(unitNames(i));
    xlim([1 T]);
    box on;
    grid on;
    ax.FontName = "Times New Roman";
    ax.LineWidth = 1.5;
    ax.FontSize = font;

    % Head
    ax = nexttile((i-1)*2 + 2);
    p_band = fill([x fliplr(x)],[head_p10 fliplr(head_p90)], ...
        band_color,"EdgeColor","none","FaceAlpha",0.45); hold on;
    plot(x,head_diu,"Color",diu_color,"LineWidth",3.5);
    plot(x,head_ddu,"Color",ddu_color,"LineWidth",3.5);
    p_ref = plot(x,head_mean,"k--","LineWidth",1.2);

    if i == 1
        title("Head (m)");
        h_band_leg = p_band;
        h_ref_leg = p_ref;
    end
    if i == n_units
        xlabel("Time (hr)");
    end

    xlim([1 T]);
    ylim([0 sp.max_h]);
    box on;
    grid on;
    ax.FontName = "Times New Roman";
    ax.LineWidth = 1.5;
    ax.FontSize = font;
end


lgd = legend([h_diu_leg h_ddu_leg h_band_leg h_ref_leg], ...
    ["DIU mean","DDU mean","SOC p10-p90","SOC mean"], ...
    "Orientation","horizontal","Box","off","FontSize",font-2);

lgd.Layout.Tile = "south";


% save figure at high resolution
fig = gcf;
fig.PaperPositionMode = 'auto';
outputFile = fullfile(pwd, "vol_highres.png");
print(fig, outputFile, "-dpng", "-r300");  % 300 DPI
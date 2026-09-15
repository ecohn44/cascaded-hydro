clear; clc; close all;

inflowFiles = [
    "/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/bonneville.csv"
    "/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/dalles.csv"
    "/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/johnday.csv"
    "/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/mcnary.csv"
];

names = ["BON","TDA","JDA","MCN"];
inflowCols = names + "_Flow_In_Inst__6Hours_0_RFC_FCST_kcfs_";

figure;
tiledlayout(2, 1); 

% Top subplot
ax1 = nexttile;
colsPlot = lines(4);
hLinesTop = gobjects(4,1);

for k = 1:4
    D = readtable(inflowFiles(k));
    t = datetime(D.("DateTime"));
    q = D.(inflowCols(k)) * 28.316846592;

    keep = year(t) >= 2018 & year(t) <= 2025 & ( ...
        (month(t) == 9 & day(t) >= 7) | ...
        month(t) == 10 | month(t) == 11 | ...
        (month(t) == 12 & day(t) <= 5));

    h = hours(t(keep) - datetime(year(t(keep)), 9, 7)) + 1;
    H = table(h, q(keep), VariableNames=["h","q"]);
    H = H(~isnan(H.q), :);

    S = groupsummary(H, "h", {@mean, @(z) prctile(z,10), @(z) prctile(z,90)}, "q");

    fill([S.h; flipud(S.h)], [S.fun2_q; flipud(S.fun3_q)], colsPlot(k,:), ...
        "FaceAlpha", 0.15, "EdgeColor", "none"); 
    hold on;

    hLinesTop(k) = plot(S.h, S.fun1_q, "Color", colsPlot(k,:), "LineWidth", 1.8);
end

ylabel("Inflow (m^3/s)", fontSize=16);
grid off;
box off;
xlim([1 2160]);

xticks([1 360 576 936 1320 1680 2160]);
xticklabels(["Sep 7", "Sep 22", "Oct 1", "Oct 16", "Nov 1", "Nov 16", "Dec 5"]);
xtickangle(-30);

%legend(ax1, hLinesTop, names, "Location", "southoutside", "Box", "off", "NumColumns", 4, "FontSize", 12);
% Bottom subplot
ax2 = nexttile;

files = [
    "/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/bonneville-soc.csv"
    "/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/dalles-soc.csv"
    "/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/johnday-soc.csv"
    "/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/mcnary-soc.csv"
];

cols = [
    "BON_Elev_Forebay_Inst_1Hour_0_CBT_REV_ft_"
    "TDA_Elev_Forebay_Inst_1Hour_0_CBT_REV_ft_"
    "JDA_Elev_Forebay_Inst_1Hour_0_CBT_REV_ft_"
    "MCN_Elev_Forebay_Inst_1Hour_0_CBT_REV_ft_"
];

names = ["BON", "TDA", "JDA", "MCN"];

soc = table();

for i = 1:4
    T = readtable(files(i));
    t = datetime(T.("DateTime"));

    keep = year(t) >= 2018 & year(t) <= 2025 & ( ...
        (month(t) == 9 & day(t) >= 7) | ...
        month(t) == 10 | month(t) == 11 | ...
        (month(t) == 12 & day(t) <= 5));

    Ti = table();
    Ti.datetime = t(keep);
    Ti.(names(i)) = T.(cols(i))(keep) * 0.3048;  % ft to m

    if i == 1
        soc = Ti;
    else
        soc = innerjoin(soc, Ti, "Keys", "datetime");
    end
end
yr = year(soc.datetime);

h = hours(soc.datetime - datetime(yr, 9, 7)) + 1;

T = table(h, soc.BON, soc.TDA, soc.JDA, soc.MCN, VariableNames=["h", "BON", "TDA", "JDA", "MCN"]);

% compute mean, 10th and 90th percentiles by hour for each reservoir
funHandles = {@mean, @(z) prctile(z,10), @(z) prctile(z,90)};
S = groupsummary(T, "h", funHandles);

% minimum forebay levels
bottom = [70, 155, 257, 330] * 0.3048;  % ft to m

% plot mean and shaded p10-p90 for each reservoir
colsPlot = lines(4);
namesPlot = ["BON","TDA","JDA","MCN"];
hLines = gobjects(4,1);
hPatches = gobjects(4,1);
for k = 1:4
    meanVec = S.("fun1_" + namesPlot(k)) - bottom(k);
    p10Vec  = S.("fun2_" + namesPlot(k)) - bottom(k);
    p90Vec  = S.("fun3_" + namesPlot(k)) - bottom(k);
    xh = S.h;
    % shaded band (create patch object)
    hPatches(k) = fill([xh; flipud(xh)], [p10Vec; flipud(p90Vec)], colsPlot(k,:), ...
        "FaceAlpha", 0.15, "EdgeColor", "none"); hold on;
    % mean line
    hLines(k) = plot(xh, meanVec, "Color", colsPlot(k,:), "LineWidth", 1.8);
end

ylabel("Head (m)", fontSize=16);
grid off;
box off;

xlim([1 2160]);
xticks([1 360 576 936 1320 1680 2160]);
xticklabels(["Sep 7", "Sep 22", "Oct 1", "Oct 16", "Nov 1", "Nov 16", "Dec 5"]);
xtickangle(-30);

% Create combined legend entries: one for shaded region and one for mean line per reservoir
legendEntries = strings(4,1);
legendHandles = gobjects(4,1);
for k = 1:4
    legendHandles(k) = hLines(k);
    legendEntries(k) = namesPlot(k);
end
% move legend for bottom subplot directly next to it (outside its axes)
lgBottom = legend(ax2, legendHandles, legendEntries, "Location", "southoutside", "Box", "off", "NumColumns", 4, "FontSize",12);

% save figure at high resolution
fig = gcf;
fig.PaperPositionMode = 'auto';
outputFile = fullfile(pwd, "figure_highres.png");
print(fig, outputFile, "-dpng", "-r300");  % 300 DPI
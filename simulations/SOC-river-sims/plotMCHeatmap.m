clear; clc; close all;

kappa_case = "range";   % "standard" or "range"

switch kappa_case
    case "standard"
        load("resultsBonferroni/monteCarloResultsk1MC.mat","results")

        years = results.years;
        thetas = results.thetas;
        P = results.p;

        m_diu = find(results.frameworks == "diu",1);
        m_ddu = find(results.frameworks == "ddu",1);
        k = find(results.kappa == 1,1);

        Y = numel(years);
        H = numel(thetas);

        power_increase = nan(Y,H);
        valid_count = zeros(Y,H);

        for y = 1:Y
            for h = 1:H
                P_diu = P(:,:,y,h,m_diu,k,:);
                P_ddu = P(:,:,y,h,m_ddu,k,:);

                valid = squeeze(all(~isnan(P_diu),[1 2]) & ...
                                all(~isnan(P_ddu),[1 2]));

                E_diu = squeeze(sum(P_diu,[1 2]));
                E_ddu = squeeze(sum(P_ddu,[1 2]));

                power_increase(y,h) = 100*(mean(E_ddu(valid))-mean(E_diu(valid))) ...
                                           / mean(E_diu(valid));
                valid_count(y,h) = sum(valid);
            end
        end

        [~,order] = sort(results.mean_inflow,"ascend");

        Z = power_increase(order,:);
        valid_count = valid_count(order,:);
        xvals = thetas;
        yvals = 1:Y;
        ylabels = string(years(order));

        xlabel_text = "\theta";
        ylabel_text = "Historical Year: Driest to Wettest";
        title_text = "DDU Generation Increase over DIU, \kappa = 1 (%)";

    case "range"
        load("resultsBonferroni/monteCarloResultskrangeMC.mat","results");

        thetas = results.thetas;
        kappa = results.kappa;
        P = results.p;

        m_diu = find(results.frameworks == "diu",1);
        m_ddu = find(results.frameworks == "ddu",1);

        H = numel(thetas);
        K = numel(kappa);

        power_increase = nan(H,K);
        valid_count = zeros(H,K);

        for h = 1:H
            for k = 1:K
                P_diu = P(:,:,:,h,m_diu,k,:);
                P_ddu = P(:,:,:,h,m_ddu,k,:);

                valid = squeeze(all(~isnan(P_diu),[1 2]) & ...
                                all(~isnan(P_ddu),[1 2]));

                E_diu = squeeze(sum(P_diu,[1 2]));
                E_ddu = squeeze(sum(P_ddu,[1 2]));

                power_increase(h,k) = 100*(mean(E_ddu(valid))-mean(E_diu(valid))) ...
                                           / mean(E_diu(valid));
                valid_count(h,k) = sum(valid,"all");
            end
        end

        Z = power_increase;
        xvals = kappa;
        yvals = thetas;
        ylabels = string(thetas);

        xlabel_text = "Forecast Error Multiplier, \kappa";
        ylabel_text = "\theta";
        title_text = "Mean DDU Generation Increase over DIU (%)";
end

figure("Color","w");

imagesc(xvals,yvals,Z);
set(gca,"YDir","normal");

lo = min(Z,[],"all","omitnan");
hi = max(Z,[],"all","omitnan");

if lo == hi
    lo = lo - 1;
    hi = hi + 1;
end

clim([lo hi]);
colormap(redblue(256));

cb = colorbar;
cb.Label.String = "Generation Increase (%)";
cb.Label.FontName = "Times New Roman";

ax = gca;
ax.XTick = xvals;
ax.YTick = yvals;
ax.YTickLabel = ylabels;
ax.FontName = "Times New Roman";
ax.FontSize = 12;
ax.LineWidth = 1.2;

xlabel(xlabel_text);
ylabel(ylabel_text);
title(title_text);

for r = 1:size(Z,1)
    for c = 1:size(Z,2)
        if isnan(Z(r,c)), continue; end

        zscaled = (Z(r,c) - lo) / (hi - lo);

        if zscaled < 0.25 || zscaled > 0.75
            txt_color = "w";
        else
            txt_color = "k";
        end

        text(xvals(c), yvals(r), sprintf("%+.2f",Z(r,c)), ...
            "HorizontalAlignment","center", ...
            "VerticalAlignment","middle", ...
            "Color",txt_color, ...
            "FontName","Times New Roman", ...
            "FontSize",12);
    end
end

box on;

disp("Number of fully converged DIU-DDU pairs in each heatmap cell:");
disp(valid_count);

function cmap = redblue(n)
    positions = [0 0.5 1];

    colors = [0.85 0.25 0.25;   % red: lowest value
              1.00 1.00 1.00;   % white: middle value
              0.05 0.25 0.65];  % blue: highest value

    cmap = interp1(positions,colors,linspace(0,1,n));
end
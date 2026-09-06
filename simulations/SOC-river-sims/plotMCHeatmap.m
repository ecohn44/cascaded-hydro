clear; clc;

load("resultsBonferroni/monteCarloResultsk1.mat","results");

years = results.years;
thetas = results.thetas;  

P = results.p;

m_diu = find(results.frameworks == "diu");
m_ddu = find(results.frameworks == "ddu");
k     = find(results.kappa == 1,1);

Y = length(years);
H = length(thetas);

power_increase = zeros(Y,H);

for y = 1:Y
    for h = 1:H

        % Dimensions retained here:
        % unit × time × year × theta × framework × kappa × scenario
        P_diu = P(:,:,y,h,m_diu,k,:);
        P_ddu = P(:,:,y,h,m_ddu,k,:);

        % Total generation for each Monte Carlo scenario
        E_diu_scen = squeeze(sum(P_diu,[1 2]));
        E_ddu_scen = squeeze(sum(P_ddu,[1 2]));

        % Average total generation across all scenarios
        E_diu = mean(E_diu_scen);
        E_ddu = mean(E_ddu_scen);

        power_increase(y,h) = 100*(E_ddu - E_diu)/E_diu;
    end
end

% Rank years from driest to wettest
[~,order] = sort(results.mean_inflow,"ascend");

power_increase = power_increase(order,:);
years_ranked   = years(order);

% Consolidated heatmap
figure("Color","w");
imagesc(thetas,1:Y,power_increase);

ax = gca;
ax.YDir = "normal";
ax.XTick = thetas;
ax.YTick = 1:Y;
ax.YTickLabel = years_ranked;
ax.FontName = "Times New Roman";
ax.FontSize = 11;

xlabel("\theta");
ylabel("Historical Year: Driest to Wettest");
title("DDU Generation Increase over DIU, \kappa = 1 (%)");

cb = colorbar;
cb.Label.String = "Generation Increase (%)";

colormap(turbo);
box off;

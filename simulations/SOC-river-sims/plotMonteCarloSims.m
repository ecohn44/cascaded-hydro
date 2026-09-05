clear; clc; close all;

load("resultsBonferroni/monteCarloResults.mat","results");

years = 2018:2025;  
kappa  = results.kappa;
failed = results.failed;
P      = results.p;

h = 1;                         % theta = 0
Y = length(years);
K = length(kappa);
power_increase = nan(Y,K);

for y = 1:Y
    for k = 1:K
        P_diu = P(:,:,y,h,1,k,:);
        P_ddu = P(:,:,y,h,2,k,:);

        E_diu = mean(squeeze(sum(P_diu,[1 2])));
        E_ddu = mean(squeeze(sum(P_ddu,[1 2])));

        power_increase(y,k) = 100*(E_ddu-E_diu)/E_diu;
    end
end

% Order years from driest to wettest
[~,order] = sort(results.mean_inflow);
power_increase = power_increase(order,:);
years = years(order);

% Heatmap
figure;
imagesc(kappa,1:Y,power_increase);
set(gca,"YTick",1:Y,"YTickLabel",years);
xlabel("Forecast Error Multiplier, \kappa");
ylabel("Historical Year: Driest to Wettest");
title("DDU Generation Increase over DIU (%)");
colorbar;

failure_rate = 100*squeeze(mean(results.failed(:,1,:,:,:),5));

figure;
tiledlayout(1,2);

for m = 1:2
    nexttile;
    plot(years,squeeze(failure_rate(:,m,:)),"-o");
    title(upper(results.frameworks(m)));
    xlabel("Year");
    ylabel("Failure Rate (%)");
    xticks(years);
    ylim([0 100]);
    grid on;
end

legend("\kappa = " + string(results.kappa),"Location","eastoutside");
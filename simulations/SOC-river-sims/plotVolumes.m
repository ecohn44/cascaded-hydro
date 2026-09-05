clear; clc; close all;

load("resultsBonferroni/monteCarloResults.mat","results");

results.years = 2018:2025;

k = find(results.kappa == 1);
h = 1;                              % theta = 0
T = size(results.V,2);
n_units = size(results.V,1);

m_diu = find(results.frameworks == "diu");
m_ddu = find(results.frameworks == "ddu");

V_min = [results.sysparams.min_V];
V_max = [results.sysparams.max_V];

for y = 1:length(results.years)

    V_diu = squeeze(mean(results.V(:,:,y,h,m_diu,k,:),7,"omitnan"));
    V_ddu = squeeze(mean(results.V(:,:,y,h,m_ddu,k,:),7,"omitnan"));

    std_diu = squeeze(mean(results.std(:,:,y,h,m_diu,k,:),7,"omitnan"));
    std_ddu = squeeze(mean(results.std(:,:,y,h,m_ddu,k,:),7,"omitnan"));

    figure;
    tiledlayout(2,2);

    for i = 1:n_units
        nexttile;

        yyaxis left
        hold on;
        a = plot(1:T,V_diu(i,:),"b","LineWidth",1.5);
        b = plot(1:T,V_ddu(i,:),"Color",[0.85 0.33 0.10],"LineWidth",1.5);
        yline(V_min(i),"r--");
        yline(V_max(i),"r--");
        ylabel("Volume");

        yyaxis right
        c = plot(1:T,std_diu(i,:),"b:","LineWidth",1.2);
        d = plot(1:T,std_ddu(i,:),":","Color",[0.85 0.33 0.10],"LineWidth",1.2);
        ylabel("Forecast Standard Deviation");

        title("Unit " + i);
        xlabel("Hour");
        grid on;

        if i == 1
            legend([a b c d], ...
                "DIU volume","DDU volume","DIU std","DDU std");
        end
    end

    sgtitle("Volume and Forecast Uncertainty: " + results.years(y));
end
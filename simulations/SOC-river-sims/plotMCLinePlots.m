clear; clc; close all;

load("resultsBonferroni/monteCarloResultsk1MC.mat","results");

P = results.p;
thetas = results.thetas;

m_diu = find(results.frameworks == "diu",1);
m_ddu = find(results.frameworks == "ddu",1);
k = find(results.kappa == 1,1);

unitNames = ["McNary","John Day","The Dalles","Bonneville"];

[~,order] = sort(results.mean_inflow,"ascend");

yearGroups = {order(1:2), order(end-1:end)};
groupNames = ["Driest years", "Wettest years"];

figure("Color","w");
tiledlayout(2,4,"TileSpacing","compact","Padding","compact");

for g = 1:2
    years_use = yearGroups{g};

    for i = 1:4
        nexttile; hold on; box on;

        for m = [m_diu m_ddu]

            mu = nan(size(thetas));
            lo = nan(size(thetas));
            hi = nan(size(thetas));

            for h = 1:numel(thetas)

                A = P(:,:,:,h,m_diu,k,:);
                B = P(:,:,:,h,m_ddu,k,:);

                valid = squeeze(all(~isnan(A),[1 2]) & all(~isnan(B),[1 2]));

                use = false(size(valid));
                use(years_use,:) = valid(years_use,:);

                E = squeeze(sum(P(i,:,:,h,m,k,:),2));   % year x scenario
                E = E(use);

                mu(h) = mean(E,"omitnan");
                lo(h) = min(E,[],"omitnan");
                hi(h) = max(E,[],"omitnan");
            end

            if m == m_diu
                c = [0 0.4470 0.7410];
                name = "DIU";
            else
                c = [0.8500 0.3250 0.0980];
                name = "DDU";
            end

            fill([thetas fliplr(thetas)], [lo fliplr(hi)], c, ...
                "EdgeColor","none", ...
                "FaceAlpha",0.15, ...
                "HandleVisibility","off");

            plot(thetas,mu,"-o", ...
                "Color",c, ...
                "MarkerFaceColor",c, ...
                "LineWidth",2, ...
                "DisplayName",name);
        end

        title(unitNames(i));
        xlabel("\theta");

        if i == 1
            ylabel(groupNames(g) + newline + "Total Power Generation");
        end

        set(gca, ...
            "FontName","Times New Roman", ...
            "FontSize",11, ...
            "XTick",thetas);
        xtickangle(0); % ensure x-axis labels are not rotated
    end
end

legend("Location","best");
sgtitle("Total Generation by Tracking Penalty, \kappa = 1");
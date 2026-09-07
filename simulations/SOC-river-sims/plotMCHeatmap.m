clear; clc; close all;

kappa_case = "standard";   % "standard" or "range"

switch kappa_case
    case "standard"
        % load("resultsBonferroni/monteCarloResultsk1.mat","results");
        load("resultsBonferroni/monteCarloResultsMC.mat","results")

        years = results.years;
        thetas = results.thetas;
        P = results.p;

        m_diu = find(results.frameworks == "diu");
        m_ddu = find(results.frameworks == "ddu");
        k = find(results.kappa == 1,1);

        Y = length(years);
        H = length(thetas);
        power_increase = zeros(Y,H);
        valid_count = zeros(Y,H);

        for y = 1:Y
            for h = 1:H
                P_diu = P(:,:,y,h,m_diu,k,:);
                P_ddu = P(:,:,y,h,m_ddu,k,:);

                % Keep only scenarios where both methods converged for all
                % units and every time step
                valid = squeeze(all(~isnan(P_diu),[1 2]) & ...
                                all(~isnan(P_ddu),[1 2]));

                E_diu_scen = squeeze(sum(P_diu,[1 2]));
                E_ddu_scen = squeeze(sum(P_ddu,[1 2]));

                E_diu = mean(E_diu_scen(valid));
                E_ddu = mean(E_ddu_scen(valid));

                power_increase(y,h) = 100*(E_ddu-E_diu)/E_diu;
                valid_count(y,h) = sum(valid);
            end
        end

        % Rank years from driest to wettest
        [~,order] = sort(results.mean_inflow,"ascend");
        power_increase = power_increase(order,:);
        valid_count = valid_count(order,:);
        years_ranked = years(order);

        figure("Color","w");
        im = imagesc(thetas,1:Y,power_increase);
        im.AlphaData = ~isnan(power_increase);

        ax = gca;
        ax.YDir = "normal";
        ax.XTick = thetas;
        ax.YTick = 1:Y;
        ax.YTickLabel = years_ranked;
        ax.FontName = "Times New Roman";
        ax.FontSize = 11;
        ax.Color = [0.75 0.75 0.75];

        xlabel("\theta");
        ylabel("Historical Year: Driest to Wettest");
        title("DDU Generation Increase over DIU, \kappa = 1 (%)");

    case "range"
        load("resultsBonferroni/monteCarloResultskrange.mat","results");

        thetas = results.thetas;
        kappa = results.kappa;
        P = results.p;

        m_diu = find(results.frameworks == "diu");
        m_ddu = find(results.frameworks == "ddu");

        H = length(thetas);
        K = length(kappa);
        power_increase = zeros(H,K);
        valid_count = zeros(H,K);

        for h = 1:H
            for k = 1:K
                P_diu = P(:,:,:,h,m_diu,k,:);
                P_ddu = P(:,:,:,h,m_ddu,k,:);

                % Keep only year-scenario pairs where both methods converged
                % for all units and every time step
                valid = squeeze(all(~isnan(P_diu),[1 2]) & ...
                                all(~isnan(P_ddu),[1 2]));

                E_diu_all = squeeze(sum(P_diu,[1 2]));
                E_ddu_all = squeeze(sum(P_ddu,[1 2]));

                E_diu = mean(E_diu_all(valid));
                E_ddu = mean(E_ddu_all(valid));

                power_increase(h,k) = 100*(E_ddu-E_diu)/E_diu;
                valid_count(h,k) = sum(valid,"all");
            end
        end

        figure("Color","w");
        im = imagesc(kappa,thetas,power_increase);
        im.AlphaData = ~isnan(power_increase);

        ax = gca;
        ax.YDir = "normal";
        ax.XTick = kappa;
        ax.YTick = thetas;
        ax.FontName = "Times New Roman";
        ax.FontSize = 11;
        ax.Color = [0.75 0.75 0.75];

        xlabel("Forecast Error Multiplier, \kappa");
        ylabel("\theta");
        title("Mean DDU Generation Increase over DIU (%)");
end

cb = colorbar;
cb.Label.String = "Generation Increase (%)";

colormap(flipud(summer));
box off;

disp("Number of fully converged DIU-DDU pairs in each heatmap cell:");
disp(valid_count);

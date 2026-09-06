clear; clc; close all;


load("resultsBonferroni/monteCarloResultsk1.mat","results");

results.years = 2018:2024;
thetas = [0, 5, 10, 15];  
test_year = 2024;

k = find(results.kappa == 1);
m_diu = find(results.frameworks == "diu");
m_ddu = find(results.frameworks == "ddu");

T = size(results.V,2);
n_units = size(results.V,1);
Y = length(results.years);
x = 1:T;
V_min = [results.sysparams.min_V]';
V_max = [results.sysparams.max_V]';


% Select driest year, wettest year, and 2025 test year
train_idx = find(results.years ~= test_year);
[~,order] = sort(results.mean_inflow(train_idx));
plot_idx = [train_idx(order(1)), train_idx(order(end)), find(results.years == test_year)];
% row_names = ["Driest", "Wettest", "Test"];

% Plot colors
diu_color = [0.00 0.45 0.74];
ddu_color = [0.85 0.33 0.10];
band_color = [0.80 0.80 0.80];


for i = 1:n_units
    figure("Color","w","Position",[100 100 1450 720]);
    tl = tiledlayout(Y,length(thetas),"TileSpacing","compact", ...
        "Padding","compact");

    V_p10  = results.SOC_p10(i,:);
    V_p90  = results.SOC_p90(i,:);
    V_mean = results.SOC_mean(i,:,7);


    for y = 1:Y

        for h = 1:length(thetas)
            ax = nexttile;
            hold on;

            V_diu = reshape(mean(results.V(i,:,y,h,m_diu,k,:),7,"omitnan"),1,T);
            V_ddu = reshape(mean(results.V(i,:,y,h,m_ddu,k,:),7,"omitnan"),1,T);

            p_band = fill([x fliplr(x)],[V_p10 fliplr(V_p90)],band_color, "EdgeColor","none","FaceAlpha",0.45);
            p_diu = plot(x,V_diu,"Color",diu_color,"LineWidth",1.5);
            p_ddu = plot(x,V_ddu,"Color",ddu_color,"LineWidth",1.5);
            p_ref = plot(x,V_mean,"k--","LineWidth",1.2);

            xlim([1 T]);
            ylim([0 results.sysparams(i).max_V]);
            grid on;
            box off;
            ax.FontName = "Times New Roman";
            ax.FontSize = 10;
            ax.GridAlpha = 0.15;

            if y == 1
                title("\theta = " + thetas(h));
            end
            if h == 1
                ylabel("Volume" + results.years(y));
            end
            if y == Y
                xlabel("Hour");
            else
                ax.XTickLabel = [];
            end

            if y == 1 && h == 1
                legend_ax = ax;
                legend_handles = [p_diu p_ddu p_band p_ref];
            end
        end
    end

    title(tl,"Unit " + i,"FontWeight","bold");
    lgd = legend(legend_ax,legend_handles, ...
        ["DIU mean","DDU mean","SOC p10-p90","SOC mean"], ...
        "Orientation","horizontal","Box","off");
    lgd.Layout.Tile = "south";
end



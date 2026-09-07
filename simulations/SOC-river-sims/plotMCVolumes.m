clear; clc; close all;

%% Select plotting case
plot_case = "subset";     % "all" or "subset"
kappa_value = 1;

results_folder = "resultsBonferroni";
all_file   = fullfile(results_folder,"monteCarloResultsk1.mat");
train_file = fullfile(results_folder,"monteCarloResultsValidation_no_umin.mat");
test_file  = fullfile(results_folder,"monteCarloResultsTest_no_umin.mat");

switch plot_case
    case "all"
        S = load(all_file,"results");
        results = S.results;

        years = results.years;
        plot_idx = 1:length(years);
        plot_results = repmat({results},1,length(years));
        row_names = "Year: " + string(years);
        thetas = results.thetas;

    case "subset"
        S = load(train_file,"results");
        train_results = S.results;

        S = load(test_file,"results");
        test_results = S.results;

        % Rank only the training years by mean streamflow
        [~,order] = sort(train_results.mean_inflow);
        wet_idx = order(end);
        dry_idx = order(4);

        years = [train_results.years(wet_idx), ...
                 train_results.years(dry_idx), ...
                 test_results.years(1)];

        plot_idx = [wet_idx,dry_idx,1];
        plot_results = {train_results,train_results,test_results};
        row_names = ["Wet train year: " + years(1), ...
                     "Dry train year: " + years(2), ...
                     "Test year: " + years(3)];

        % Plot theta values available in both files
        thetas = intersect(train_results.thetas, ...
                           test_results.thetas,"stable");
end

R = length(plot_idx);
T = size(plot_results{1}.V,2);
n_units = size(plot_results{1}.V,1);
x = 1:T;

% Plot colors
diu_color = [0.00 0.45 0.74];
ddu_color = [0.85 0.33 0.10];
band_color = [0.80 0.80 0.80];

for i = 1:n_units
    figure("Color","w","Position",[100 100 1450 720]);
    tl = tiledlayout(R,length(thetas),"TileSpacing","compact", ...
        "Padding","compact");

    for r = 1:R
        results = plot_results{r};
        y = plot_idx(r);

        k = find(results.kappa == kappa_value,1);
        m_diu = find(results.frameworks == "diu",1);
        m_ddu = find(results.frameworks == "ddu",1);

        V_mean = reshape(results.SOC_mean(i,:,y),1,T);
        V_p10  = reshape(results.SOC_p10(i,:,y),1,T);
        V_p90  = reshape(results.SOC_p90(i,:,y),1,T);

        for h = 1:length(thetas)
            theta_idx = find(results.thetas == thetas(h),1);

            ax = nexttile;
            hold on;

            V_diu = reshape(mean(results.V(i,:,y,theta_idx,m_diu,k,:),7),1,T);
            V_ddu = reshape(mean(results.V(i,:,y,theta_idx,m_ddu,k,:),7),1,T);

            p_band = fill([x fliplr(x)],[V_p10 fliplr(V_p90)], ...
                band_color,"EdgeColor","none","FaceAlpha",0.45);
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

            if r == 1
                title("\theta = " + thetas(h));
            end

            if h == 1
                ylabel({row_names(r),"Volume"});
            end

            if r == R
                xlabel("Hour");
            else
                ax.XTickLabel = [];
            end

            if r == 1 && h == 1
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

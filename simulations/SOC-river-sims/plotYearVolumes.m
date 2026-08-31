function plotYearVolumes(results_folder)
% plotYearVolumes("resultsOracle")

close all;
plot_years = false;

years = 2018:2025; 
thetas = [0.01 0.1 1];
T = 2160;

for a = 1:length(thetas)

    theta = thetas(a);

    % Load each year
    for k = 1:length(years)
        folder = fullfile(results_folder, ...
            sprintf('%d_dry_det_T%d',years(k),T));

        file = fullfile(folder, ...
            sprintf('results_theta%d.mat',round(100*theta)));

        S = load(file,'X','sysparams');
        V(:,:,k) = S.X(:,1:5:end)';
    end

    n = size(V,1);
    T = size(V,2);
    colors = lines(length(years));

    figure('Position',[100 100 1200 650]);
    tiledlayout(2,2,'TileSpacing','compact');

    for i = 1:n
        nexttile;
        hold on;

        Y = squeeze(V(i,:,:));

        if plot_years
            % Plot each year
            for k = 1:length(years)
                h(k) = plot(Y(:,k), 'Color',colors(k,:), 'LineWidth',1.2);
            end

            labels = arrayfun(@num2str,years, ...
                'UniformOutput',false);

        else
            % Plot only min-max range and mean
            Vmin  = min(Y,[],2);
            Vmax  = max(Y,[],2);
            Vmean = mean(Y,2);

            h_band = fill( ...
                [1:T fliplr(1:T)], ...
                [Vmin' fliplr(Vmax')], ...
                [0.75 0.85 1], ...
                'FaceAlpha',0.4, ...
                'EdgeColor','none');

            h_mean = plot(Vmean,'k--','LineWidth',1.5); hold on
            h_min = plot(Vmin,'b','LineWidth',1);
            h_max = plot(Vmax,'b','LineWidth',1);

            h = [h_band h_mean h_min h_max];
            labels = {'Min-max range','Mean', 'Min', 'Max'};
        end

        title(sprintf('Unit %d: %s', i, S.sysparams(i).name));
        xlabel('Time Step');
        ylabel('Volume');
        ylim([0 S.sysparams(i).max_V]);
        grid on;
    end

    lgd = legend(h,labels);
    lgd.Layout.Tile = 'east';

    sgtitle(sprintf('Volume Trajectories: \\theta = %g',theta));

    clear V h
end

end
function plotThetaVolumes(results_folder)

files = dir(fullfile(results_folder,'results_theta*.mat'));

% Load results
for k = 1:length(files)
    load(fullfile(files(k).folder,files(k).name), 'X','theta', 'sysparams');
    V(:,:,k) = X(:,1:5:end)';   % unit x time x theta
    thetas(k) = theta;
end

% Sort theta values
[thetas,idx] = sort(thetas);
V = V(:,:,idx);

n = size(V,1);
T = size(V,2);
colors = turbo(length(thetas));

figure('Position',[100 100 1200 650]);
tiledlayout(2,2,'TileSpacing','compact');

for i = 1:n
    nexttile;
    hold on;

    Y = squeeze(V(i,:,:));       % time x theta
    Vmin  = min(Y,[],2);
    Vmax  = max(Y,[],2);
    Vmean = mean(Y,2);

    % Min-max shaded region
    h_band = fill( ...
        [1:T, fliplr(1:T)], ...
        [Vmin', fliplr(Vmax')], ...
        [0.75 0.85 1], ...
        'FaceAlpha',0.4, ...
        'EdgeColor','none');

    % Individual theta trajectories
    for k = 1:length(thetas)
        h_theta(k) = plot(Y(:,k), 'Color',colors(k,:), 'LineWidth',1.2);
    end

    % Mean trajectory
    h_mean = plot(Vmean,'k--','LineWidth',1);
    title(sprintf('Unit %d: %s', i, sysparams(i).name));
    ylim([0 sysparams(i).max_V]);
    xlim([0, T])
    xlabel('Time Step');
    ylabel('Volume');
    grid on;
end

% Single global legend
labels = [{'Min-max range'}, ...
    arrayfun(@(x) sprintf('\\theta = %g',x), ...
    thetas,'UniformOutput',false), ...
    {'Mean'}];

lgd = legend([h_band,h_theta,h_mean],labels);
lgd.Layout.Tile = 'east';

sgtitle('Volume Tracking-Penalty Sensitivity');

end
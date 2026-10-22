function plotRiskSpread(alphas_MC)

    [T, n, K] = size(alphas_MC);
    tt = 1:T;

    % Mean colors
    cMean{1} = [0.0000, 0.4470, 0.7410];   % blue
    cMean{2} = [0.8500, 0.3250, 0.0980];   % orange
    cMean{3} = [0.9290, 0.6940, 0.1250];   % yellow

    % Lighter fill colors for shadow bands
    cFill{1} = [0.65, 0.82, 0.97];
    cFill{2} = [0.98, 0.75, 0.60];
    cFill{3} = [0.98, 0.90, 0.60];

    figure('Color','w');
    hold on; box on; grid on;

    h = gobjects(n+1,1);

    % --- shaded IQR band + mean line for each unit ---
    for i = 1:n
        A = squeeze(alphas_MC(:,i,:));   % T x K

        alpha_mean = median(A, 2, 'omitnan');
        alpha_q25  = prctile(A', 25)';   % T x 1
        alpha_q75  = prctile(A', 75)';

        % Shadow band
        fill([tt, fliplr(tt)], ...
             [alpha_q25', fliplr(alpha_q75')], ...
             cFill{i}, ...
             'FaceAlpha', 0.35, ...
             'EdgeColor', 'none', ...
             'HandleVisibility', 'off');

        % Mean line
        h(i) = plot(tt, alpha_mean, ...
            'Color', cMean{i}, ...
            'LineWidth', 2.8);
    end

    % --- Bonferroni reference line ---
    h(n+1) = yline(1/n, '--', ...
        'Color', [0.3 0.3 0.3], ...
        'LineWidth', 2.0);

    xlabel('Time Step', 'FontSize', 16);
    ylabel('\alpha_i', 'FontSize', 16);

    legend(h, {'Unit 1 Mean', 'Unit 2 Mean', 'Unit 3 Mean', 'Bonferroni: 1/n'}, ...
        'Location', 'best', ...
        'FontSize', 13, ...
        'Box', 'on');

    xlim([1 T]);
    ylim([0 0.5]);   % tighter axis to emphasize convergence
    set(gca, ...
        'FontSize', 15, ...
        'LineWidth', 1.2, ...
        'Box', 'on');

end
function plotRiskSpread(alphas_MC)

    [T, n, K] = size(alphas_MC);
    tt = 1:T;

    % Mean colors
    cMean{1} = [0.0000, 0.4470, 0.7410];   % dark blue
    cMean{2} = [0.8500, 0.3250, 0.0980];   % orange/red
    cMean{3} = [0.9290, 0.6940, 0.1250];   % yellow

    % Light trajectory colors
    cLite{1} = [0.65, 0.82, 0.97];         % light blue
    cLite{2} = [0.98, 0.75, 0.60];         % light orange
    cLite{3} = [0.98, 0.90, 0.60];         % light yellow

    figure('Color','w'); 
    hold on; box on; grid on;

    % --- plot MC trajectories first ---
    for i = 1:n
        for k = 1:K
            plot(tt, alphas_MC(:,i,k), ...
                'Color', cLite{i}, ...
                'LineWidth', 0.75, ...
                'HandleVisibility','off');
        end
    end

    % --- plot mean trajectories on top ---
    h = gobjects(n,1);
    for i = 1:n
        alpha_mean = mean(alphas_MC(:,i,:), 3, 'omitnan');
        h(i) = plot(tt, alpha_mean, ...
            'Color', cMean{i}, ...
            'LineWidth', 2.2);
    end

    xlabel('Time Step');
    ylabel('\alpha_i');
    title('Monte Carlo Risk Allocation Trajectories');
    legend(h, {'Unit 1 Mean','Unit 2 Mean','Unit 3 Mean'}, ...
        'Location','best');

    xlim([1 T]);
    ylim([0 1]);   % remove if needed
    set(gca,'FontSize',14);

end
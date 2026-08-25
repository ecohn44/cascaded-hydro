function plotStreamflows(q)
    % q: T+lag x n matrix, each column is a streamflow time series

    [T, n] = size(q);
    t = (1:T)';        

    figure;
    for i = 1:n
        subplot(n, 1, i);
        plot(t, q(:, i), 'LineWidth', 3);
        ylim([0 1.1*max(q(:,i))]); 
        xlim([1, T]);
        
        ylabel(sprintf('q_%d', i), 'FontSize', 16);
        set(gca, 'FontSize', 16); 

        if i == 1
            title('Streamflow Time Series', 'FontSize', 20);
        end
        if i == n
            xlabel('Time (hour)', 'FontSize', 16);
        else
            set(gca, 'XTickLabel', []);  % hide x labels for middle plots
        end

        grid on;
    end
end

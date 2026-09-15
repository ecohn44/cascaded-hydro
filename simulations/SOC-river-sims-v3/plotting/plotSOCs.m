function plotSOCs(soc)

    [T, n] = size(soc);
    t = (1:T)';        

    figure;
    for i = 1:n
        subplot(n, 1, i);
        plot(t, soc(:, i), 'LineWidth', 3);
        ylim([0 1.1*max(soc(:,i))]); 
        xlim([1, T]);
        
        ylabel(sprintf('V_%d', i), 'FontSize', 16);
        set(gca, 'FontSize', 16); 

        if i == 1
            title('Historical Volume Time Series', 'FontSize', 20);
        end
        if i == n
            xlabel('Time (hour)', 'FontSize', 16);
        else
            set(gca, 'XTickLabel', []);  % hide x labels for middle plots
        end

        grid on;
    end
end

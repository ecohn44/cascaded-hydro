function plotSSH(phi_vals, alpha_vals, eps)

    T = size(alpha_vals, 1);
    tt = 1:1:T;
    n = size(alpha_vals, 2);

    % Create a new figure
    f = figure('Name','SSH Diagnostics','NumberTitle','off');

    % Manually adjust phi(t = 1) since under init conditions 
    phi_vals(1) = 1;
    
    % Figure 1: Pr(V in [Vmin, Vmax]) 
    subplot(3, 1, 1, 'Parent', f); 
    plot(tt, phi_vals, 'LineWidth', 1.5, 'DisplayName', 'System Reliability');
    yline(1-eps, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Reliability Bound');
    title('Total System Reliability');
    ylabel('\phi'); 
    xlabel('Time (hr)');
    xlim([1, T]);
    ylim([0, 1.1]);
    legend('Location', 'west');
    grid on;
    
    % Figure 2: Risk Attribution Weights 
    subplot(3, 1, 2, 'Parent', f);  
    plot(tt, alpha_vals, 'LineWidth', 1.5);
    title('Risk Attribution Weights');
    xlabel('Time (hr)');
    ylabel('\alpha_i'); 
    xlim([1, T])
    grid on;
    
    % Create a cell array of strings for the legend labels
    legendLabels = cell(1, n);
    for i = 1:n
        legendLabels{i} = sprintf('Unit %d', i);
    end
    
    % Apply the legend
    legend(legendLabels, 'Location', 'west');

    % Figure 3: Effective Risk Budget
    subplot(3, 1, 3, 'Parent', f);  
    plot(tt, alpha_vals*eps, 'LineWidth', 1.5);
    yline(eps/n, 'r--', 'LineWidth', 1.5);
    title('Effective Risk Budget');
    xlabel('Time (hr)');
    ylabel('\epsilon_i'); 
    xlim([1, T]);

    % Create a cell array of strings for the legend labels
    legendLabels = cell(1, n+1);
    for i = 1:n
        legendLabels{i} = sprintf('\\epsilon^{SSH}_{%d}', i);
    end
    legendLabels{n+1} = '\epsilon^{BON}';
    
    legend(legendLabels, 'Location', 'west');
    grid on;


end
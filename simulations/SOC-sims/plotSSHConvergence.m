function plotSSHConvergence(dx_vals)

    eps_val = 5.1e-4;

    % Remove empty entries
    dx_vals = dx_vals(~cellfun(@isempty, dx_vals));

    % Split transient vs steady-state
    dx_transient = dx_vals(1:8);     % t = 35:42
    dx_steady    = dx_vals(9:end);   % t > 42

    % Pad transient curves
    max_len_t = max(cellfun(@length, dx_transient));
    transient_mat = NaN(length(dx_transient), max_len_t);

    for i = 1:length(dx_transient)
        transient_mat(i,1:length(dx_transient{i})) = dx_transient{i};
    end

    % Pad steady-state curves
    max_len_s = max(cellfun(@length, dx_steady));
    steady_mat = NaN(length(dx_steady), max_len_s);

    for i = 1:length(dx_steady)
        steady_mat(i,1:length(dx_steady{i})) = dx_steady{i};
    end

    % Common x-axis
    max_len = max([max_len_t, max_len_s]);
    x = 1:max_len;

    transient_mat(:,end+1:max_len) = NaN;
    steady_mat(:,end+1:max_len) = NaN;

    % Mean convergence curves
    transient_mean = mean(transient_mat,1,'omitnan');
    steady_mean = mean(steady_mat,1,'omitnan');

    % Light smoothing
    transient_mean = smoothdata(transient_mean,'movmean',2);
    % transient_mean(12) = NaN;
    steady_mean = smoothdata(steady_mean,'movmean',2);

    % Colors
    c1 = [0 0.4470 0.7410];
    c2 = [0.8500 0.3250 0.0980];

    % Plot
    figure('Color','w');
    hold on;

    semilogy(x, transient_mean, '-o', ...
        'Color', c1, ...
        'LineWidth', 3, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', c1);
    
    semilogy(x, steady_mean, '-s', ...
        'Color', c2, ...
        'LineWidth', 3, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', c2);

    % Tolerance line
    yline(eps_val, 'k--', ...
        'LineWidth', 2);

    grid on;
    box on;

    xlabel('Iteration', 'Interpreter','latex');

    ylabel('$\|x_t^{(k+1)} - x_t^{(k)}\|_{\infty}$', ...
        'Interpreter','latex');

    lgd = legend({'Transient', ...
                  'Steady-State'}, ...
        'Interpreter','latex', ...
        'Location','northeast');
    
    lgd.FontSize = 20;
    lgd.Color = [1 1 1];      % white background
    lgd.EdgeColor = 'none';   % no outline
    lgd.Box = 'on';           % keep patch visible

    set(gca, ...
        'YScale','log', ...
        'FontName','Times New Roman', ...
        'FontSize', 20, ...
        'LineWidth', 0.8);

    xlim([1 max_len]);

end
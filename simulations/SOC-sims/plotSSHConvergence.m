function plotSSHConvergence(dx_vals, eps_val)

    if nargin < 2
        eps_val = 1e-3;
    end

    % Remove empty entries
    dx_nonempty = dx_vals(~cellfun(@isempty, dx_vals));

    % Pad with NaNs
    max_len = max(cellfun(@numel, dx_nonempty));
    dx_mat = NaN(numel(dx_nonempty), max_len);

    for i = 1:numel(dx_nonempty)
        v = dx_nonempty{i};
        dx_mat(i,1:numel(v)) = v;
    end

    % Mean convergence
    dx_mean = mean(dx_mat, 1, 'omitnan');

    % Plot
    figure('Color','w');
    semilogy(1:max_len, dx_mean, '-o', ...
        'LineWidth', 2.5, ...
        'MarkerSize', 6);

    hold on;
    yline(eps_val, '--', ...
        'LineWidth', 2);

    grid on;
    box on;

    xlabel('Iteration', 'Interpreter','latex');
    ylabel('$\|x^{(k+1)} - x^{(k)}\|_{\infty}$', 'Interpreter','latex');

    set(gca, ...
        'FontName','Times New Roman', ...
        'FontSize', 20, ...
        'LineWidth', 0.8);

    xlim([1 max_len]);
end
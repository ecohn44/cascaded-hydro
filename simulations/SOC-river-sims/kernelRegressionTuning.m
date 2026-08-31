function [best_W, best_sigma, mean_rmse, fold_rmse] = kernelRegressionTuning(inflow_norm, soc_ref, W_grid, sigma_grid, years)
% Tune W and sigma by leave-one-scenario-out cross-validation.


    % Extract problem dimensions
    [n, T, S] = size(inflow_norm); % [4, 2160, 8]
    fold_rmse = zeros(numel(W_grid), numel(sigma_grid), S);
    avg_rmse = zeros(numel(W_grid), numel(sigma_grid), S);

    % For each window size 
    for i = 1:numel(W_grid)
        % For each smoothness parameter
        for j = 1:numel(sigma_grid)
            for test = 1:S

                % Get training year indices
                train = setdiff(1:S, test);

                % Get historical inflow data 
                inflow_test = inflow_norm(:, :, test);

                % Extract test SOC reference 
                soc_hat = zeros(n, T);

                % Initialize with inital volume state 
                soc_hat(:, 1) = soc_ref(:, 1, test);

                % Generate reference for t = 2:T
                for t = 1:T-1

                    % Calculate next best SOC state from training data 
                    soc_hat(:, t + 1) = kernelRegression( ...
                        permute(inflow_norm(:, :, train), [2 1 3]), ...
                        permute(soc_ref(:, :, train), [2 1 3]), ...
                        inflow_test(:, 1:t)', ...
                        soc_hat(:, 1:t)', ...
                        W_grid(i), sigma_grid(j));
                end

                % Average of the six training trajectories
                soc_average = mean(soc_ref(:, :, train), 3);
                error_avg = soc_average(:, 2:end) - soc_ref(:, 2:end, test);
                avg_rmse(i, j, test) = sqrt(mean(error_avg(:).^2));

                % Calcuate test error 
                error = soc_hat(:, 2:end) - soc_ref(:, 2:end, test);
                fold_rmse(i, j, test) = sqrt(mean(error(:).^2));

                % Plot simulated and reference SOCs
                plotSOCs(soc_hat, soc_ref(:, :, test), soc_average, years(test));
            end
        end
    end

    % Average RMSE over all years to get mean RMSE for each parameter pair
    mean_rmse = mean(fold_rmse, 3);

    % Find index of minimum RMSE (converted to column vector)
    [~, k] = min(mean_rmse(:));

    % Convert linear index to subscripts 
    [i, j] = ind2sub(size(mean_rmse), k);

    % Identify optimal parameters
    best_W = W_grid(i);
    best_sigma = sigma_grid(j);

    fprintf('Best W = %g, best sigma = %.4g, mean RMSE = %.4f\n', ...
        best_W, best_sigma, mean_rmse(i, j));

    fprintf('Average SOC RMSE: %.4f\n', mean(average_rmse));

end


function plotSOCs(soc_hat, soc_s, soc_avg, year)

    [n, T] = size(soc_hat);
    t = (1:T)';        

    figure;
    for i = 1:n
        subplot(n, 1, i);
        plot(t, soc_hat(i, :), 'LineWidth', 3); hold on 
        plot(t, soc_s(i, :), 'r--', 'LineWidth', 2); hold on 
        plot(t, soc_avg(i, :), 'g--', 'LineWidth', 2);
        xlim([1, T]);
        
        ylabel(sprintf('V_%d', i), 'FontSize', 16);
        set(gca, 'FontSize', 16); 

        if i == 1
            title('SOC Trajectories ' + string(year), 'FontSize', 20);
        end
        if i == n
            xlabel('Time (hour)', 'FontSize', 16);
        else
            set(gca, 'XTickLabel', []);  % hide x labels for middle plots
        end

        grid on;
    end
    
end

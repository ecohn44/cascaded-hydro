function LMP = simulatePrice(T, n, enable_inflation)
%   Inputs:
%       T    - Total duration of the simulation in hours.
%       n    - Number of units (columns in output).
%       enable_inflation - Boolean to trigger price spike.
%
%   Outputs:
%       LMP  - T x n matrix of price signals.
%              (Normalized per unit: Min=1, Max=1.5)
    
    % Default enable_inflation to false
    if nargin < 3
        enable_inflation = false;
    end

    seed = 42;
    % Set Random Seed
    rng(seed);
    % Generate Base 24-Hour Profile
    % Anchors: Start(0), Low(4am), Peak(8am), SolarLow(1pm), Peak(5pm), End(24)
    t_anchors = [0,    4,   8,    13,   17,   24];
    p_anchors = [0.5, 0.4, 1.0,  0.2,  1.0,  0.5];
    
    hours_in_day = 1:24;
    base_day = interp1(t_anchors, p_anchors, hours_in_day, 'pchip');
    
    % Repeat for T Hours to get Base Vector
    num_days = ceil(T / 24);
    base_full = repmat(base_day, 1, num_days);
    base_signal = base_full(1:T)'; % Transpose to make it a column vector (Tx1)
    
    % Create Matrix for n Units and Add Noise
    % We replicate the base signal into n columns
    LMP_matrix_raw = repmat(base_signal, 1, n);
    
    % Generate noise (Uniform distribution centered on 0)
    % Adjustment factor: 0.1 means noise is roughly +/- 5% of magnitude
    noise_level = 0.15; 
    noise = (rand(T, n) - 0.5) * noise_level;
    
    % Add noise to the base signal
    LMP_noisy = LMP_matrix_raw + noise;
    
    % Normalize Each Unit Independently
    target_min = 1.0;
    target_max = 1.5;
    
    LMP = zeros(T, n); % Pre-allocate
    
    for i = 1:n
        col = LMP_noisy(:, i);
        min_val = min(col);
        max_val = max(col);
        
        % Normalize column i
        LMP(:, i) = target_min + (col - min_val) .* ...
                    ((target_max - target_min) / (max_val - min_val));
    end

    % Apply Inflation to First Two Units
    if enable_inflation
        start_idx = max(1, T - 11);
        LMP(start_idx:end, 1:min(n, 2)) = LMP(start_idx:end, 1:min(n, 2)) * 2;
    end
    
    t = (1:T)';
    % Plotting (Stacked Subplots)
    figure('Name', 'LMP Simulation per Unit', 'Color', 'w');
    
    for i = 1:n
        subplot(n, 1, i);
        plot(t, LMP(:, i), 'LineWidth', 1.5);
        grid on;
        
        % Styling
        ylabel('Price ($/MWh)');
        title(['Unit ', num2str(i), ' LMP']);
        xlim([1 T]);
        ylim([0.9 max(LMP(:)) + 0.1]); % Adjusted to fit inflation
        
        % Add peak markers for visual check
        xline(8:24:T, ':b', 'Alpha', 0.3); % 8am markers
        xline(17:24:T, ':r', 'Alpha', 0.3); % 5pm markers
        
        % Only put x-label on the bottom plot
        if i == n
            xlabel('Time (Hours)');
        else
            set(gca, 'XTickLabel', []); % Hide numbers for cleaner stack
        end
    end
end
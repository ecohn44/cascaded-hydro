function [soc_next, weights] = kernelRegression(inflow_norm, soc_ref, inflow_history, soc_history, window, sigma)

% Historical data:
%   inflow_norm  T x n x S   normalized inflow for S scenarios
%   soc_ref      T x n x S   oracle SOC; 
%
% Data available online at time t:
%   inflow_history  t x n
%   soc_history     t x n     measured SOC, including the current SOC
%
% Output:
%   soc_next        1 x n     reference for time t+1
%   weights         S x 1     similarity weight for each scenario

    [T, n, S] = size(inflow_norm);
    t = size(inflow_history, 1);

    assert(t < T, 'No t+1 reference exists when t >= T.');

    % Create look back window 
    idx = max(1, t - window + 1):t;

    % Real time look back window 
    x = [reshape(inflow_history(idx, :), [], 1); reshape(soc_history(idx, :), [], 1)];

    similarity = zeros(S, 1);
    % For each historical scenario s in S
    for s = 1:S
        % Create historical look back window 
        y = [reshape(inflow_norm(idx, :, s), [], 1); reshape(soc_ref(idx, :, s), [], 1)];
        % Calculate the squared distance for similarity measure
        distance2 = sum((x - y).^2);
        % Calculate the similarity for the current scenario
        similarity(s) = exp(-distance2 / (2 * sigma^2 * numel(x)));
    end

    total_similarity = sum(similarity);
    if total_similarity == 0
        weights = ones(S, 1) / S;
    else
        weights = similarity / total_similarity;
    end
    candidates = reshape(soc_ref(t + 1, :, :), n, S)';
    soc_next = weights' * candidates;
end

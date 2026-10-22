function R = estimateR(T, n, lag, q, m)

    alpha0 = m.AR_const;
    alpha1 = m.AR_coef;

    % Preallocate forecast and residuals
    q_hat = nan(T, n);
    E     = nan(T, n);
    
    % AR(1) forecast (q already starts at lag)
    for t = 1:T
        % Forecast streamflow based on lagged observation 
        q_hat(t,:) = alpha0 + alpha1*q(t,:);

        % Calculate forecast error based on actual observation 
        E(t,:)     = q(t+lag,:) - q_hat(t,:);
    end

    % Compute correlation of forecast errors
    R = corr(E, 'Rows', 'complete');

end
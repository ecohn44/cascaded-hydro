function std_hat = forecast_error(t, kappa, q_error, up_release, framework, model, sys)

    n = numel(sys);
    std_hat = model.AR_std * ones(1, n);

    if framework == "det"
        std_hat = zeros(1, n);
        return
    end

    for i = 1:n
        if framework == "ddu" && i > 1 && t > 1
            var_hat = model.omega + model.alpha*(q_error(i)^2) + model.gamma*up_release(i);
            std_hat(i) = kappa * sqrt(var_hat);
        end
    end

end

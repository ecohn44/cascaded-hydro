function std_hat = forecast_error(t, kappa, q_error, up_release, up_ramp, framework, model, sys)

    n = numel(sys);
    std_hat = zeros(n,1);

    if framework == "det"
        return
    end

    % McNary has constant forecast uncertainty
    std_hat(1) = 9.0; % kfs

    for i = 2:n
        j = i - 1;  % model 1=JDA, 2=TDA, 3=BON

        if framework == "diu" || t == 1
            std_hat(i) = model(j).AR_std;

        else % DDU ARCH-X model
            release_scaled = max(0,(up_release(i-1) - model(j).release_min)/model(j).release_range);
            ramp_scaled = max(0,(abs(up_ramp(i-1)) - model(j).ramp_min)/model(j).ramp_range);

            var_hat = model(j).omega + model(j).alpha*q_error(i)^2 + model(j).gamma_release*release_scaled + model(j).gamma_ramp*ramp_scaled;
            std_hat(i) = kappa*sqrt(max(var_hat,0));
        end
    end
end
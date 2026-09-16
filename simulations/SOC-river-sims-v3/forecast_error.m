function std_hat = forecast_error(t, kappa, q_error, up_release, framework, model, sys)

    n = numel(sys);
    std_hat = zeros(n,1);

    if framework == "det"
        return
    end

    % McNary has constant forecast uncertainty
    std_hat(1) = kappa*model(1).AR_std;

    for i = 2:n
        j = i - 1;  % model 1=JDA, 2=TDA, 3=BON

        if framework == "diu" || t == 1
            std_hat(i) = kappa*model(j).AR_std;
        
        else % ddu garch model 
            var_hat = model(j).omega + model(j).alpha*q_error(i)^2 + model(j).gamma*up_release(i-1);
            std_hat(i) = kappa*sqrt(max(var_hat,0));
        end
    end

end
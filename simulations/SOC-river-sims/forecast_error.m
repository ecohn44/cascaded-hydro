% Helper: Construct inflow forecast error 
function std_hat = forecast_error(t, q_error, up_release, framework, model, sys)
    
    n = numel(sys);
    std_hat = model.AR_std*ones(1,n);    % estimated forecast variance 

    % Calculate inflow for each unit
    for i = 1:n
       
        if framework == "det"
           % No error estimation
           std_hat =  zeros(1,n);

        elseif framework == "ddu"
            if i > 1 && t > 1
                % Forecast conditional variance using GARCH-X
                var_hat_norm =  model.omega + model.alpha*(q_error(i)^2) + model.gamma*(up_release(i)); 
                std_hat(i) = sqrt(var_hat_norm);
            end
        end
    end
end
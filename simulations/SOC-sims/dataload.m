%% ========================================================================
% DATALOAD FUNCTIONS
% ========================================================================

function [params, sysparams] = dataload(N)
    % Outputs:
    %   params  - Fitted OLS model parameters for forecasting inflow
    %   sysparams - Reservoir operation parameters
    

    %% Seasonal Forecasting Parameters
    params(1) = struct( ...
        'season', 'dry', ...    % Hydrology season
        'constant', -0.001, ... % DDU Params
        'coef1', 0.83, ...      % DDU inflow_lag1 [m3/hr]
        'coef2', 0.17, ...      % DDU outflow_lag [m3/hr]
        'AR_const', 0.0020, ... % DIU params
        'AR_coef', 0.950, ...   % DIU inflow_lag1
        'AR_std', 0.00375, ...  % Updated to be 5% of flow (0.075)
        'omega', 0.00001, ... % GARCH-X params
        'alpha', 0.008, ... % GARCH-X params
        'gamma', 0.0004);   % GARCH-X params

    % Copy for Wet Season
    params(2) = params(1);
    params(2).season = 'wet';
    params(2).gamma = 0.0008;

    % Unit 01
    sysparams(1) = struct( ...
        'unit', 1, ...
        'a', 5.0, ...  % Hydraulic head coef (max effective voltage)
        'b', 0.45, ... % Hydraulic head coef (concavity)
        'min_ut', .05, ...
        'max_ut', 0.10, ...
        'min_h', 0, ...
        'max_h', 5.0, ...
        'min_V', 0, ...
        'max_V', 1, ...
        'RR_dn', -0.25, ...
        'RR_up', 0.0025, ...
        'F', 1, ...
        'SOC', 0.5);

    % Copy for Unit 02
    sysparams(2) = sysparams(1); 
    sysparams(2).unit = 2;
    sysparams(2).SOC = 0.05;


    %% Compute Piecewise Linear Approximation
    for i = 1:numel(sysparams)
        % Compute derived parameters for this unit
        d = sysparams(i);
        
        % Create sub-intervals and midpoints
        [lb, rb, ref] = create_intervals_and_references(d.min_h, d.max_h, N);
        
        % Store in the struct
        sysparams(i).h_lbounds = lb;
        sysparams(i).h_rbounds = rb;
        sysparams(i).h_refvals = ref;

        % Calculate initial volume 
        sysparams(i).V0 = d.SOC*(d.max_V - d.min_V) + d.min_V;
    end
end

function m3s = cfs_to_m3s(cfs)
    % Conversion factor: 1 cfs is approximately 0.0283168 m³/s
    conversion_factor = 0.0283168;
    m3s = cfs * conversion_factor;
end


% Function to create N sub-intervals and midpoint references for PWL
function [left_bounds, right_bounds, reference_values] = create_intervals_and_references(min_h, max_h, N)
    % Create N+1 breakpoints to define N intervals
    breakpoints = linspace(min_h, max_h, N+1);  % vector of breakpoints

    % Create interval bounds
    left_bounds  = breakpoints(1:end-1);  % first N points
    right_bounds = breakpoints(2:end);    % last N points

    % Create reference vector with midpoints
    reference_values = (left_bounds + right_bounds) / 2;
end
%% ========================================================================
% DATALOAD FUNCTIONS
% ========================================================================

function [params, sysparams, scenarioparams] = dataload(n, N)
    % Inputs:
    %   n        - number of hydropower units
    %   N        - number of PWL segments for head approximation
    %
    % Outputs:
    %   params    - Fitted OLS / GARCH model parameters for forecasting inflow
    %   sysparams - Reservoir operation parameters (1 x n struct array)
    

    %% Seasonal Forecasting Parameters
    params(1) = struct( ...
        'season',   'dry', ...    % Hydrology season
        'constant', 0.0020, ...   % DDU Params 
        'coef1',    0.950,  ...   % DDU inflow_lag1 [m3/hr] 
        'coef2',    0,  ...       % DDU outflow_lag [m3/hr] 
        'AR_const', 0.0020, ...   % DIU params
        'AR_coef',  0.950, ...    % DIU inflow_lag1
        'AR_std',   0.006, ...   % DIU residual variance % 0.0028
        'omega',    (0.006^2), ... % GARCH-X params 
        'alpha',    0.008,   ...      % GARCH-X params % 0.8
        'gamma',    0.003);        % GARCH-X params %6.6e-4

    params(2) = struct( ...
        'season',   'wet', ...    % Hydrology season
        'constant', 0.00, ...     % DDU Params
        'coef1',    0.995,  ...   % DDU inflow_lag1 [m3/hr]
        'coef2',    0,  ...       % DDU outflow_lag [m3/hr] 
        'AR_const', 0.00, ...     % DIU params
        'AR_coef',  0.995, ...    % DIU inflow_lag1
        'AR_std',   0.0035, ...   % DIU residual variance 
        'omega',    (0.0035^2), ... % GARCH-X params 
        'alpha',    0.8,   ...      % GARCH-X params 
        'gamma',    0.001);        % GARCH-X params 


    %% Base Unit Template 
    base_unit = struct( ...
        'unit',   1, ...      
        'name',   'Unit 01', ... 
        'a',      5.0, ...    % Hydraulic head coef (max effective voltage)
        'b',      0.45, ...   % Hydraulic head coef (concavity)
        'min_ut', 0.01, ...   
        'max_ut', 0.05, ...
        'min_h',  0, ...
        'max_h',  5.0, ...
        'min_V',  0, ...
        'max_V',  1, ...
        'RR_dn',  -0.015, ...
        'RR_up',  0.01, ...
        'F',      1, ...
        'SOC',    0.5);       % Initial start of charge


    %% Build sysparams for n units
    sysparams = repmat(base_unit, 1, n);

    for i = 1:n
        % Numeric unit ID
        sysparams(i).unit = i;

        % String name
        sysparams(i).name = sprintf('Unit %02d', i); 
    end


    %% Compute Piecewise Linear Approximation and initial volume V0
    for i = 1:numel(sysparams)
        % Compute derived parameters for this unit
        d = sysparams(i);
        
        % Create sub-intervals and midpoints
        [lb, rb, ref] = create_intervals_and_references(d.min_h, d.max_h, N);
        
        % Store in the struct
        sysparams(i).h_lbounds = lb;
        sysparams(i).h_rbounds = rb;
        sysparams(i).h_refvals = ref;

        % Calculate initial volume from SOC
        sysparams(i).V0 = d.SOC*(d.max_V - d.min_V) + d.min_V;
    end

    %% Drought scenario parameters
    scenarioTemplate = struct( ...
        'mode',         '',  ...
        'q0',           [], ...  
        'amp1',         [],  ...
        'amp2',         [],  ...
        'w1',           [],  ...
        'w2',           [],  ...
        't0',           [],  ...
        'nEvents',      [],  ...
        'daysPerEvent', [],  ...
        'tauHours',     [],  ...
        'unitDelay',    [],  ...
        'startSteps',   [],  ...
        'recoverHours', []);

    % Preallocate 1x2 struct array with same fields
    scenarioparams(1) = scenarioTemplate;
    scenarioparams(2) = scenarioTemplate;
    scenarioparams(3) = scenarioTemplate;

    % Pulse-type event (uses t0 as FRACTIONS of horizon)
    scenarioparams(1).mode         = 'pulse';
    scenarioparams(1).amp1         = 0.4;      % 40% drop in inflow (dry) or +40% (wet)
    scenarioparams(1).amp2         = 0.3;      % 30% drop / bump for second pulse
    scenarioparams(1).w1           = 8;        % first drought pulse lasts 8 hours
    scenarioparams(1).w2           = 4;        % second drought pulse lasts 4 hours
    scenarioparams(1).t0           = [0.3, 0.5];  % pulses at 30% and 80% of horizon
    scenarioparams(1).unitDelay    = 12;       % time units between drought events 
    scenarioparams(1).startSteps   = 2;        % time before initial drought event begins 

    % Extended event (single-decay events)
    scenarioparams(2).mode         = 'extended';
    scenarioparams(2).amp1         = 0.4;    % Magnitude of drought event
    scenarioparams(2).nEvents      = 1;      % number of events
    scenarioparams(2).daysPerEvent = 1.2;    % drought length 
    scenarioparams(2).tauHours     = 12;     % decay rate
    scenarioparams(2).unitDelay    = 16;     % time units between drought events 
    scenarioparams(2).startSteps   = 2;      % time before initial drought event begins 
    scenarioparams(2).recoverHours = 10;     % recovery time from drought event to q0

    scenarioparams(3).mode         = 'constant';

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

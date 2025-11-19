%% ========================================================================
% DATALOAD FUNCTIONS
% ========================================================================

function [params, sysparams, droughtparams] = dataload(n, N)
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
        'constant', -0.001, ...   % DDU Params
        'coef1',    0.83,  ...    % DDU inflow_lag1 [m3/hr]
        'coef2',    0.17,  ...    % DDU outflow_lag [m3/hr]
        'AR_const', 0.0020, ...   % DIU params
        'AR_coef',  0.950, ...    % DIU inflow_lag1
        'AR_std',   0.00375, ...  % Updated to be 5% of flow (0.075)
        'omega',    0.00001, ...  % GARCH-X params
        'alpha',    0.008,   ...  % GARCH-X params
        'gamma',    0.0004);      % GARCH-X params

    % Copy for Wet Season
    params(2) = params(1);
    params(2).season = 'wet';
    params(2).gamma  = 0.0008;


    %% Base unit template (from original Unit 01)
    base_unit = struct( ...
        'unit',   1, ...      % will be overwritten per unit
        'name',   'Unit 01', ... % will be overwritten per unit
        'a',      5.0, ...    % Hydraulic head coef (max effective voltage)
        'b',      0.45, ...   % Hydraulic head coef (concavity)
        'min_ut', 0.05, ...
        'max_ut', 0.10, ...
        'min_h',  0, ...
        'max_h',  5.0, ...
        'min_V',  0, ...
        'max_V',  1, ...
        'RR_dn',  -0.25, ...
        'RR_up',  0.0025, ...
        'F',      1, ...
        'SOC',    0.5);       % will be set per unit (0.5 or 0.05)


    %% Build sysparams for n units
    sysparams = repmat(base_unit, 1, n);

    for i = 1:n
        % Numeric unit ID
        sysparams(i).unit = i;

        % String name
        sysparams(i).name = sprintf('Unit %02d', i);

        % SoC rule:
        %   Units 1..n-1 -> SOC = 0.5
        %   Unit n       -> SOC = 0.05
        if i < n
            sysparams(i).SOC = 0.5;
        else
            sysparams(i).SOC = 0.45;
        end
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
    droughtTemplate = struct( ...
        'mode',         '',  ...
        'amp1',         [],  ...
        'amp2',         [],  ...
        'w1',           [],  ...
        'w2',           [],  ...
        't0',           [],  ...
        'nEvents',      [],  ...
        'daysPerEvent', [],  ...
        'tauHours',     [] );

    % Preallocate 1x2 struct array with same fields
    droughtparams(1) = droughtTemplate;
    droughtparams(2) = droughtTemplate;
    droughtparams(3) = droughtTemplate;

    % Pulse-type drought (uses t0 as FRACTIONS of horizon; T+lag handled later)
    droughtparams(1).mode         = 'pulse';
    droughtparams(1).amp1         = 0.4;      % 40% drop in inflow (dry) or +40% (wet)
    droughtparams(1).amp2         = 0.3;      % 30% drop / bump for second pulse
    droughtparams(1).w1           = 8;        % first drought pulse lasts 8 hours
    droughtparams(1).w2           = 4;        % second drought pulse lasts 4 hours
    droughtparams(1).t0           = [0.3, 0.8];  % pulses at 30% and 80% of horizon
    % Fields not used by pulse mode remain [] (nEvents, daysPerEvent, tauHours)

    % Extended drought (single-decay events)
    droughtparams(2).mode         = 'extended';
    droughtparams(2).amp1         = 0.2;
    droughtparams(2).nEvents      = 3;
    droughtparams(2).daysPerEvent = 2;
    droughtparams(2).tauHours     = 36;

    droughtparams(3).mode         = 'constant';

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

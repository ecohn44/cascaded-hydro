function q = droughtSimulator(T, lag, season, mode, droughtParams)
% droughtSimulator  Generate inflow under different drought scenarios.
%
% INPUTS
%   q0            baseline inflow (scalar or vector length T+lag)
%   T, lag        horizon and lag
%   season        'wet' → positive pulses, 'dry' → negative pulses
%   mode          'pulse'    -> localized drought (or flood) pulses
%                 'extended' -> extended drought regime
%   droughtParams struct containing parameters for the chosen mode
%
%
% OUTPUT
%   q  : (T+lag) × 1 inflow time series (nonnegative)
%
% -------------------------------------------------------------------------

    n = T + lag;
    q0 = droughtParams.q0;

    % Build baseline q
    q = q0 * ones(n,1);


    % Dispatch based on mode
    switch lower(mode)
        case 'constant'
            % Do nothing: q stays equal to the baseline q0
        case 'pulse'
            q = applyPulseDrought(q, T, lag, season, droughtParams);

        case 'extended'
            q = applyExtendedDrought(q, T, lag, season, droughtParams);
    end

    % Enforce nonnegativity
    q = max(q, 0);

end


function q = applyPulseDrought(q, T, lag, season, droughtParams)

    n = T + lag;

    % Unpack all params (we assume they exist)
    t0  = droughtParams.t0(:)';  
    amp1 = droughtParams.amp1;
    amp2 = droughtParams.amp2;
    w1   = droughtParams.w1;
    w2   = droughtParams.w2;

    % Map fractional or absolute t0 → valid indices
    idx = @(x) round((x < 1) .* (lag + x*T) + (x >= 1) .* x);
    t_idx = max(1, min(n, arrayfun(idx, t0)));

    % Season modifier
    if strcmpi(season, 'dry')
        sgn = -1;
    else
        sgn = 1;
    end

    % Loop over pulses
    for k = 1:numel(t_idx)
        if k == 1
            amp = sgn * amp1;
            wid = round(w1);
        else
            amp = sgn * amp2;
            wid = round(w2);
        end

        half = floor((wid-1)/2);
        i1   = max(1, t_idx(k) - half);
        i2   = min(n, i1 + wid - 1);

        q(i1:i2) = q(i1:i2) .* (1 + amp);
    end
end


function q = applyExtendedDrought(q, T, lag, season, droughtParams)
    % q is baseline vector (length T+lag)

    nTot = T + lag;

    % Unpack parameters (assume they are provided)
    amp1         = droughtParams.amp1;          % fractional amplitude (e.g. 0.4)
    nEvents      = droughtParams.nEvents;       % number of drought blocks
    daysPerEvent = droughtParams.daysPerEvent;  % length of each event in days

    % Optional: exponential time constant (in hours).
    % If not provided, use half the event length as a rough default.
    if isfield(droughtParams, 'tauHours')
        tauHours = droughtParams.tauHours;
    else
        tauHours = daysPerEvent * 24 / 2;  % decay mostly within event
    end

    % Convert to time steps (assume 1 step = 1 hour)
    eventLen = max(1, round(daysPerEvent * 24));  % steps per event
    tauSteps = max(1, round(tauHours));           % exponential time constant in steps

    % --- NEW: gap between events = floor(daysPerEvent/2) days ---
    % in time steps, this is roughly half the event length
    gapLen = floor(eventLen / 2);   % = floor(daysPerEvent/2 * 24)

    % Season sign: dry → negative amplitude (drop), wet → increase
    if strcmpi(season, 'dry')
        sign_mult = -1;
    else
        sign_mult =  1;
    end
    amp_eff = sign_mult * amp1;   % effective fractional change

    % Pointer to the start of the next event
    curStart = 1;

    % Loop over events
    for e = 1:nEvents
        startIdx   = curStart;
        fullEndIdx = startIdx + eventLen - 1;

        % If this event would start beyond horizon, or not fit, stop
        if startIdx > nTot || fullEndIdx > nTot
            break;
        end

        % Apply exponential drought shape over this event
        for tIdx = startIdx:fullEndIdx
            % Position within this event (0-based)
            tInEvent = tIdx - startIdx;   % 0,1,...,eventLen-1

            % Exponential decay factor:
            %   tInEvent = 0         -> factor = 1
            %   tInEvent >> tauSteps -> factor ~ 1 + amp_eff
            decayShape = 1 - exp(- double(tInEvent) / double(tauSteps));
            factor     = 1 + amp_eff * decayShape;

            % Apply to baseline q
            q(tIdx) = q(tIdx) * factor;
        end

        % Move start pointer for the next event:
        % end of this event + 1 step + gapLen steps
        curStart = fullEndIdx + gapLen + 1;
        if curStart > nTot
            break;
        end
    end
end

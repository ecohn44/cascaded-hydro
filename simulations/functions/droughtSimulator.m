function q = droughtSimulator(T, lag, season, mode, droughtParams)
% droughtSimulator  Generate inflow under different drought scenarios,
%                   with optional *spatial lag* of the drought anomaly.
%
% INPUTS
%   T, lag        : horizon and AR lag (same as before)
%   season        : 'wet' / 'dry'
%   mode          : 'constant' / 'pulse' / 'extended'
%   droughtParams : struct with at least
%                     .q0        baseline inflow (scalar)
%                   optionally:
%                     .shiftMult integer (0,1,2,...) = how many *lag*-blocks
%                                to delay the drought anomaly
%
% OUTPUT
%   q  : (T+lag) × 1 inflow time series
%
%   For example, if you set in the caller
%       dp.shiftMult = i-1;
%   then reservoir i sees the same drought shape, but delayed by
%       (i-1)*lag  time steps.

    n = T + lag;
    q0 = droughtParams.q0;

    % Baseline (no drought) is the same for all units
    base = q0 * ones(n,1);

    % Start from baseline
    q = base;

    % Shape the drought on top of baseline 
    switch lower(mode)
        case 'constant'
            % Do nothing: stays at baseline q0

        case 'pulse'
            q = applyPulseDrought(q, T, lag, season, droughtParams);

        case 'extended'
            q = applyExtendedDrought(q, T, lag, season, droughtParams);

        otherwise
            error('Unknown drought mode: %s', mode);
    end

    % Apply spatial lag to the drought
    deltaq = q - base;   % length n

    % Optional field: how many lag-blocks to delay this unit
    shiftMult = 0;
    if isfield(droughtParams, 'shiftMult')
        shiftMult = 2*droughtParams.shiftMult;
    end

    % Number of time steps to shift anomaly
    shiftSteps = round(shiftMult * lag);

    if shiftSteps > 0
        % (Downstream) fill early part with no anomaly, shift anomaly forward
        deltaq = [zeros(shiftSteps,1); deltaq(1:end-shiftSteps)];
    elseif shiftSteps < 0
        % (Upstream) drought arrives earlier
        s = -shiftSteps;
        deltaq = [deltaq(1+s:end); zeros(s,1)];
    end

    % Reconstruct inflow and enforce nonnegativity
    q = base + deltaq;
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

    % Gap between events = floor(daysPerEvent/2) days 
    gapLen = floor(eventLen);   % = floor(daysPerEvent/2)

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

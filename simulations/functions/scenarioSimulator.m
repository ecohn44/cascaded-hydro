function q = scenarioSimulator(T, lag, season, mode, scenarioParams)
% droughtSimulator  Generate inflow under different scenarios,
%                   with optional *spatial lag* of the scenario anomaly.
%
% INPUTS
%   T, lag        : horizon and AR lag (same as before)
%   season        : 'wet' / 'dry'
%   mode          : 'constant' / 'pulse' / 'extended'
%   droughtParams : struct with at least
%                     .q0        baseline inflow (scalar)
%                   optionally:
%                     .shiftMult integer (0,1,2,...) = how many *lag*-blocks
%                                to delay the scenario anomaly
%
% OUTPUT
%   q  : (T+lag) × 1 inflow time series


    n = T + lag;
    q0 = scenarioParams.q0;

    % Baseline (no drought) is the same for all units
    base = q0 * ones(n,1);

    % Start from baseline
    q = base;

    % Shape the drought on top of baseline 
    switch lower(mode)
        case 'constant'
            % Do nothing: stays at baseline q0

        case 'pulse'
            q = applyPulse(q, T, lag, season, scenarioParams);

        case 'extended'
            q = applyExtended(q, T, lag, season, scenarioParams);

        case 'ramp'
            q = applyRamp(q, T, lag, season, scenarioParams);

        otherwise
            error('Unknown drought mode: %s', mode);
    end

    % Apply spatial lag to the drought
    deltaq = q - base;   % length n

    % Number of time steps to shift drought events between units 
    if ~strcmpi(mode,'ramp')
        shiftSteps = round(scenarioParams.startSteps + scenarioParams.shiftMult * scenarioParams.unitDelay);
        deltaq     = shift_with_zeros(deltaq, shiftSteps);
    end

    % Reconstruct inflow and enforce nonnegativity
    q = base + deltaq;
    q = max(q, 0);


end


function q = applyRamp(q, T, lag, season, scenarioParams)

    % Store base flow
    q0 = q;
    
    % Create normalized ramp from 0 to +1(wet) or -1(dry)
    if season == "wet"
        endpoint = .5;
    else
        endpoint = -.5;
    end

    ramp = linspace(0,endpoint,T+lag)';   % T×1

    q = q0 + scenarioParams.k * ramp;

end


function q = applyPulse(q, T, lag, season, scenarioParams)

    n = T + lag;

    % Unpack all params (we assume they exist)
    t0  = scenarioParams.t0(:)';  
    amp1 = scenarioParams.amp1;
    amp2 = scenarioParams.amp2;
    w1   = scenarioParams.w1;
    w2   = scenarioParams.w2;

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


function q = applyExtended(q, T, lag, season, droughtParams)
    % q is baseline vector (length T+lag)

    nTot = T + lag;

    % Unpack parameters (assume they are provided)
    amp1         = droughtParams.amp1;          % fractional amplitude (e.g. 0.4)
    nEvents      = droughtParams.nEvents;       % number of drought blocks
    daysPerEvent = droughtParams.daysPerEvent;  % length of each event in days
    recLen = droughtParams.recoverHours;        % recovery length 


    % Exponential time constant 
    tauSteps= droughtParams.tauHours;

    % Time steps per event (convert to hours)
    eventLen = floor(daysPerEvent * 24);  

    % Time step gap between events 
    gapLen = floor(eventLen/4);

    % Season sign: dry → negative; wet → positive
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

            % Exponential decay factor
            decayShape = 1 - exp(- double(tInEvent) / double(tauSteps));
            factor     = 1 + amp_eff * decayShape;

            % Apply to baseline q
            q(tIdx) = q(tIdx) * factor;
        end

        % apply delayed recovery rate 
        endFactor = 1 + amp_eff * (1 - exp(- double(eventLen-1) / double(tauSteps)));
        recStart  = fullEndIdx + 1;
        recEnd    = min(nTot, fullEndIdx + recLen);

        H = recEnd - recStart;   % last tRec value (>=0)

        for tIdx = recStart:recEnd
            tRec = tIdx - recStart;  % 0...H
        
            % normalized exponential: 1 at tRec=0, 0 at tRec=H exactly
            recShape = (exp(-double(tRec)/double(tauSteps)) - exp(-double(H)/double(tauSteps))) ...
                     / (1 - exp(-double(H)/double(tauSteps)));
        
            factor = 1 + (endFactor - 1) * recShape;   % endFactor -> 1 exactly
            q(tIdx) = q(tIdx) * factor;
        end

        % Move start pointer for the next event
        curStart = fullEndIdx + gapLen + 1;
        if curStart > nTot
            break;
        end
    end
end


function y = shift_with_zeros(x, k)
% Shift vector x by k steps, padding with zeros
% k > 0 shifts forward (delays signal)
% k < 0 shifts backward (advances)

    n = numel(x);
    y = zeros(n,1);

    if k == 0
        y = x;
        return;
    end

    if abs(k) >= n
        % signal moved out of window
        return;
    end

    if k > 0
        y(k+1:end) = x(1:end-k);
    else
        s = -k;
        y(1:end-s) = x(s+1:end);
    end
end
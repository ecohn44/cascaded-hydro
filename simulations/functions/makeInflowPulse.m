function q = makeInflowPulse(q0, T, lag, t0, amp1, amp2, w1, w2, season)
    % makeInflowPulse  Generate inflow with up to two pulses
    % q0     baseline (scalar or vector length T+lag)
    % T,lag  horizon and lag
    % t0     pulse positions (fractions <1 or absolute indices)
    % amp1,amp2  fractional amplitudes (0.2 = +20%)
    % w1,w2      pulse widths (in hours)
    % season     'wet' → positive pulses, 'dry' → negative pulses
    %
    % Returns q (T+lag x 1) >= 0
    
    n = T + lag;
    
    if isscalar(q0)
        q = q0 * ones(n,1);
    else
        q = q0(:);
    end
    
    t0 = t0(:)';  % ensure row vector
    idx = @(x) round((x < 1) .* (lag + x*T) + (x >= 1) .* x);
    t_idx = max(1, min(n, arrayfun(idx, t0)));
    
    sign_mult = 1;
    if strcmpi(season, 'dry')
        sign_mult = -1;
    end
    
    for k = 1:numel(t_idx)
        if k == 1
            amp = sign_mult * max(0, amp1);
            wid = max(0, round(w1));
        else
            amp = sign_mult * max(0, amp2);
            wid = max(0, round(w2));
        end
        if amp == 0 || wid == 0
            continue
        end
        half = floor((wid-1)/2);
        i1 = max(1, t_idx(k) - half);
        i2 = min(n, i1 + wid - 1);
        q(i1:i2) = q(i1:i2) .* (1 + amp);
    end
    
    q = max(q, 0);
    end
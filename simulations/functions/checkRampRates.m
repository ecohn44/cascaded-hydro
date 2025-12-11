function checkRampRates(X, sysparams, tol)
% checkRampRates  Verify ramp-up and ramp-down constraints over a simulation.
%
%   X         : T x (5*n_units) simulation matrix
%   sysparams : struct array, one per unit, with fields .RU and .RD
%   tol       : small numerical tolerance (e.g., 1e-6)

if nargin < 3
    tol = 1e-6;
end

[T, ncols] = size(X);
n_units = ncols / 5;

fprintf('--- Ramp Rate Check ---\n');

for i = 1:n_units
    base = 5*(i-1);
    u = X(:, base + 3);    % release column for unit i

    du = diff(u);          % u(t) - u(t-1), length T-1

    RU = sysparams(i).RR_up;  % ramp-up limit per step
    RD = sysparams(i).RR_dn;  % ramp-down limit per step

    % Violations
    viol_up_idx   = find(du >  RU);      % ramp-up too big
    viol_down_idx = find(du <  RD);      % ramp-down too big

    max_up   = max(du, [], 'omitnan');
    max_down = min(du, [], 'omitnan');         % most negative change

    fprintf('Unit %d:\n', i);
    fprintf('  Max ramp up   = %+g (limit %+g)\n', max_up,  RU + tol);
    fprintf('  Max ramp down = %+g (limit -%g)\n', max_down, RD - tol);

    if isempty(viol_up_idx) && isempty(viol_down_idx)
        fprintf('  ✅ No ramp violations.\n');
    else
        fprintf('  ❌ Violations found:\n');
        if ~isempty(viol_up_idx)
            fprintf('    Ramp-up violated at t = %s\n', mat2str(viol_up_idx'+1));
        end
        if ~isempty(viol_down_idx)
            fprintf('    Ramp-down violated at t = %s\n', mat2str(viol_down_idx'+1));
        end
    end
end
end
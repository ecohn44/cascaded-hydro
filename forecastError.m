%% ========================================================================
% FORECAST ERROR MODEL
% ========================================================================

function [sg, s2] = forecastError(i, uLev, duR, s2prev, AX, omgK, kap, polK, sfloor, sceil)
% Conditional standard deviation from the ARCH-X model, one look-ahead step
%
%   sigma^2 = omega_i(kappa) + alpha_i * s2prev
%             + kappa * ( gu_i*ut + gr_i*|dut| )
%
% Inputs:
%   i       - unit index
%   uLev    - upstream release level at this step (kcfs)
%   duR     - upstream ramp at this step (kcfs/hr)
%   s2prev  - squared residual at m = 1, conditional variance for m >= 2
%   AX      - fitted ARCH-X parameters, read from data/archx_params.csv
%   omgK    - kappa-adjusted intercept, one per unit
%   kap     - decision-dependence scaling; 0 is the DIU arm
%   polK    - PolicyKappa, scales the whole returned standard deviation
%   sfloor  - lower clip on the returned standard deviation
%   sceil   - upper clip on the returned standard deviation
%
% Outputs:
%   sg      - conditional standard deviation
%   s2      - the CAPPED variance, to feed the next step's recursion
%
% s2prev follows the standard multi-step GARCH recursion: the observed
% squared residual at the first step, the previous conditional variance after
% that. Exogenous inputs are min-max normalised on the estimation range and
% clipped at zero so the model never extrapolates below it.
%
% kap = 0 gives the decision-independent arm, in which sigma is the
% unconditional standard deviation carried through the ARCH recursion only.

%% 1. Normalised exogenous inputs
ut  = max((uLev - AX.umin(i)) / AX.urng(i), 0);
dut = max((duR  - AX.rmin(i)) / AX.rrng(i), 0);

%% 2. Conditional variance
s2  = omgK(i) + AX.alpha(i)*s2prev + kap*(AX.gu(i)*ut + AX.gr(i)*dut);
s2  = max(s2, 1e-12);
sg  = sqrt(s2);

if polK ~= 1
    sg = polK*sg;  s2 = sg^2;
end

%% 3. Clip, and return the clipped variance
sg  = min(max(sg, sfloor), sceil);
s2  = sg^2;                       % feed the CAPPED variance into the recursion
end

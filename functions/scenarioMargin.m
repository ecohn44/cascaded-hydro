%% ========================================================================
% SCENARIO-PROGRAM MARGIN
% ========================================================================

function [z, info] = scenarioMargin(n, Rcorr, S, r, seed)
% The scenario-program multiplier, as an order statistic
%
%   z = scenarioMargin(n, Rcorr, S, r, seed)
%   scenarioMargin                          % no args: print the design table
%
%   WHY THIS IS A SCALAR AND NOT S EXTRA LP ROWS. The scenario approach
%   (Calafiore-Campi) samples S realisations of the shock and enforces the
%   constraint in all of them. Here the look-ahead storage row is
%
%       Vplan(i,k) + sqrt(k)*sigma_i*zeta(i,k)  >=  Vlow(i,k)
%
%   with the decision entering only through Vplan. Enforcing that for all S
%   samples is the same as enforcing it at the worst sample -- the rows
%   share no structure beyond Vplan -- so S copies collapse to one row with
%   a larger margin. And since realtimeSimulation's margin is z*sqrt(k)*sigma_i
%   while the shock is sqrt(k)*sigma_i*zeta(i,k), SIGMA AND sqrt(k) CANCEL:
%   the requirement is just z >= -zeta(i,k). The multiplier depends only on
%   (n, Rcorr, S, r) and can be computed once per run, which is why
%   'scenario' is the same one-line hook as 'robust'.
%
%   THE HORIZON DOES NOT APPEAR, AND THAT IS THE WHOLE POINT. An earlier
%   version of this file took the maximum over units AND look-ahead steps,
%   i.e. it enforced
%       P( all i, all k are simultaneously fine ) >= 1 - eps.
%   That is NOT the paper's constraint. realtimeSimulation sets nTails = n, one
%   lower tail per project, and SSH targets 1 - eps AT EVERY LOOK-AHEAD HOUR
%   separately -- a per-step joint constraint over the n projects, not a
%   joint constraint over the whole horizon. The wrong statistic returns
%   2.94 where the right one returns 2.23 at eps = 0.05, which would have
%   made the scenario benchmark a straw man: 30% more margin than Bonferroni
%   for a constraint nobody asked for. Because zeta(i,k) is standard normal
%   at every k, its distribution does not depend on k, so ONE scalar is
%   correct and is applied at every step, exactly as Bonferroni's z is.
%
%   VALIDATION. At eps = 0.05, n = 4, independent units, the exact per-step
%   multiplier is Phi^-1((1-eps)^(1/n)) = 2.2340 against Bonferroni's
%   Phi^-1(1-eps/n) = 2.2414 -- a gap of 0.0074, which is the number already
%   recorded in the project notes for the Bonferroni-to-exact gap at n = 4.
%   The discarding variant below must converge to 2.2340. If it does not,
%   something here is wrong; do not report the run.
%
%   THREE CONDITIONS make the collapse valid, all of which hold in
%   realtimeSimulation.m: one here-and-now plan with no per-scenario recourse;
%   additive uncertainty that does not multiply a decision variable; and
%   Sigma fixed within the rolling solve (paper Sec. IV, eq. 6c). Break any
%   one -- recourse especially -- and you really would need S copies. Do not
%   reuse this helper elsewhere without rechecking all three.
%
%   DISCARDING. With r = 0 the multiplier is the sample MAXIMUM, the pure
%   scenario program. It grows without bound in S because it enforces every
%   drawn sample and the Gaussian support is unbounded: more samples make it
%   MORE conservative, which is the opposite of the usual intuition and is
%   the headline finding of this benchmark. Sampling with discarding
%   (Campi-Garatti) removes the r worst samples and takes the (r+1)-th
%   largest; with r = floor(eps*S) it targets the true joint quantile.
%
%   WHAT THE DISCARDING VARIANT ACTUALLY ISOLATES, and this is why it is
%   worth running rather than just citing: it gets the joint level exactly
%   right but applies ONE standardised margin to every unit, which is equal
%   allocation. The gap between it and SSH is therefore the ALLOCATION gain
%   alone, with the union-bound conservatism already removed. That makes it
%   a sharper control than Bonferroni for the paper's central claim.
%
%   INPUTS
%     n      projects
%     Rcorr  n-by-n cross-project correlation of the hourly shock
%     S      samples
%     r      samples to discard (0 = pure scenario program)
%     seed   optional, default 20260928
%
%   OUTPUT
%     z      the multiplier to use in place of Phi^-1(1 - eps/n)
%     info   struct with the sorted statistic and the realised level

if nargin == 0, print_design_table(); return, end
if nargin < 4 || isempty(r),    r = 0;           end
if nargin < 5 || isempty(seed), seed = 20260928; end

assert(isscalar(S) && S >= 1 && S == fix(S), 'S must be a positive integer');
assert(isscalar(r) && r >= 0 && r == fix(r) && r < S, ...
    'r must be an integer in [0, S-1]; got %g with S = %d', r, S);
assert(isequal(size(Rcorr),[n n]), 'Rcorr must be %d-by-%d', n, n);

[L, p] = chol(Rcorr, 'lower');
if p > 0
    [V, D] = eig((Rcorr+Rcorr')/2);
    A = V*sqrt(diag(max(diag(D), 1e-10)));
    L = A ./ max(vecnorm(A, 2, 2), eps);       % renormalise to unit diagonal
end

rng(seed, 'twister');
w = L * randn(n, S);                           % S draws of the n-vector
M = max(-w, [], 1);                            % per-sample joint statistic
Ms = sort(M(:), 'descend');
z  = Ms(r+1);

info = struct('S',S,'r',r,'seed',seed,'sorted',Ms,'z',z, ...
    'discard_fraction',r/S,'quantile_level',1-(r+1)/(S+1), ...
    'mean_stat',mean(M),'max_stat',Ms(1));
end

% ----------------------------------------------------------------------
function print_design_table()
%PRINT_DESIGN_TABLE  What multiplier does each (S, r) choose?
%
%   Run this BEFORE committing compute. It answers the design question in
%   seconds, with no dispatch solves, and it is also the correctness check:
%   the r = eps*S column must converge to the exact value printed above it.

n = 4;  EPS = 0.05;
root = fileparts(fileparts(mfilename('fullpath')));
Rcorr = eye(n);
try
    C = dataload(root, 2024, 0.20, 7);
    if isfield(C,'Rcorr') && isequal(size(C.Rcorr),[n n])
        Rcorr = C.Rcorr;
        fprintf('using the case cross-project correlation\n');
    end
catch
    fprintf('could not load the case; using identity correlation\n');
end

zb = -sqrt(2)*erfcinv(2*(1-EPS/n));            % Bonferroni
ze = -sqrt(2)*erfcinv(2*(1-EPS)^(1/n));        % exact, independent units
fprintf('\nn = %d, eps = %.2f\n', n, EPS);
fprintf('  Bonferroni          z = %.4f\n', zb);
fprintf('  exact, independent  z = %.4f   (gap %.4f)\n', ze, zb-ze);
fprintf(['  the r = eps*S column below must converge to the exact value\n' ...
         '  (lower, if the loaded Rcorr has positive cross-correlation)\n\n']);

fprintf('  %7s %10s %8s %10s\n', 'S', 'r = 0', 'r', 'r = eps*S');
for S = [10 20 50 100 200 500 1000 5000]
    r = floor(EPS*S);
    fprintf('  %7d %10.4f %8d %10.4f\n', S, ...
        scenarioMargin(n,Rcorr,S,0), r, scenarioMargin(n,Rcorr,S,r));
end
fprintf(['\nRead the r = 0 column downward: the pure scenario program gets MORE\n' ...
         'conservative as S grows, because it enforces every drawn sample and the\n' ...
         'Gaussian support is unbounded. The r = eps*S column is the discarding\n' ...
         'variant and should settle near the exact value. Report both.\n']);
end

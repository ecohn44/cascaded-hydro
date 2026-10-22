%% ========================================================================
% CASE WINDOWING
% ========================================================================

function Cs = sliceCase(C, h0, nHours)
% Take an nHours window of a case, starting at hour h0
%
%   Cs = sliceCase(C, 4032, 168)      % the week beginning 17 Mar
%
%   realtimeSimulation has no start-hour option -- T = min([cfg.MaxHours, ...]) at
%   line 38, so every run begins at hour 1, which is 00:00 on 1 October.
%   That matters for a timing study, because October is the worst month of
%   the water year for constraint attainment.
%
%   Measured on the production full-year runs (theta 400, eps 0.05, SSH,
%   four scenarios, n = 4), the share of hours where SSH reports the joint
%   target unattainable -- no feasible point reaches 1 - eps, so the base
%   solve exits before generating a single cut:
%
%       Oct 16.1   Nov  6.2   Dec  5.9   Jan  5.6   Feb  4.3   Mar 2.8
%       Apr  2.9   May  6.3   Jun  4.6   Jul  6.7   Aug 13.3   Sep 12.0
%
%   Annual mean 7.2%. The first 168 hours -- what the smoke test was using
%   -- run at 8.6%, above average. Cleanest 168-hour windows: h4032 (17 Mar)
%   at 0.3%, h4704 (14 Apr) at 1.6%, h3360 (18 Feb) at 2.1%. Worst: h504
%   (22 Oct) at 24.7%.
%
%   So for a CPU-time study, start in March. Attainment is then a non-issue
%   at n = 4, and any rise with n is attributable to n rather than to the
%   season. Say in the caption which window was used.
%
%   CAVEAT. The first few hours of any window have clamped upstream history,
%   because C.lag reaches back up to 12 hours and realtimeSimulation clamps the
%   index at zero. That is equally true of hour 1 in every production run,
%   so a sliced window is no worse -- but do not read the first ~12 hours of
%   a short window as representative.
%
%   Slice BEFORE tiling. Row slicing and column replication commute, and
%   slicing first keeps tileCase's diagnostics on the window being run.

T = size(C.ref,1);
assert(h0 >= 1, 'h0 must be at least 1');
assert(h0 + nHours - 1 <= T, ...
    'window h0=%d + %d hours runs past the end of the case (%d hours)', ...
    h0, nHours, T);

idx = h0:(h0 + nHours - 1);
Cs  = C;

for f = {'ref','qLow','qHigh','obs','local','outObs','zseq'}
    if isfield(C,f{1}), Cs.(f{1}) = C.(f{1})(idx,:); end
end
for f = {'routedExo','upReleaseExo','lmp'}
    if isfield(C,f{1}), Cs.(f{1}) = C.(f{1})(idx); end
end

Cs.sliceStart = h0;
Cs.sliceHours = nHours;
end

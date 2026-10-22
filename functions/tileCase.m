%% ========================================================================
% SYNTHETIC CASCADE TILING
% ========================================================================

function Ct = tileCase(C, m, varargin)
% Build an m-block synthetic cascade by chaining copies of a case
%
%   Ct = tileCase(C, m)                       % default LocalMode 'matched'
%   Ct = tileCase(C, m, 'LocalMode','zeromean')
%   Ct = tileCase(C, m, 'LocalMode','raw')
%   Ct = tileCase(C, m, 'Quiet',true)
%
%   C is whatever dataload returns. The result is a single chain of
%   n = m*numel(C.units) projects:
%
%       MCN -> JDA -> TDA -> BON -> MCN_2 -> JDA_2 -> ... -> BON_m
%
%   so block b's head takes its inflow from block b-1's tail. A real chain,
%   not m independent copies.
%
%   WHAT THIS IS FOR. A synthetic scaling study for the CPU-time comment.
%   Report ONLY solve time, iteration and cut counts, and problem size.
%   Never generation, revenue or IVI as physical results.
%
%   WHY THE BLOCK IS TILED, NOT A SINGLE UNIT. SSH earns its advantage by
%   allocating risk unevenly across HETEROGENEOUS units -- at the fitted
%   WY2024 parameters The Dalles carries 52.6% of its variance in the ramp
%   term against McNary's 4.8%. Replicating one unit would make every unit
%   identical and, by Proposition 1, the allocation would collapse onto the
%   uniform split and the method would appear to do nothing.
%
%   LOCALMODE -- read this, the first version got it wrong.
%   Block 1's head is fed by routedExo, mean 124.9 kcfs, and its four units
%   add 12.7 kcfs of local inflow, so its tail releases about 138. Block 2's
%   head is fed by that 138. The 12.7 kcfs surplus is therefore inherited by
%   every downstream block and never leaves.
%     'raw'       replicates the locals untouched, so the surplus also
%                 ACCUMULATES: 12.7 kcfs more per block.
%     'zeromean'  removes the mean local in blocks 2..m. This stops the
%                 accumulation but NOT the inheritance -- block 2 still runs
%                 12.7 kcfs wetter than block 1 forever. Measured at n = 8
%                 it left block 2 sitting 9 to 21 percentage points fuller
%                 than block 1, Bonneville_2 at 82% of usable storage, with
%                 45% more spill and the base LP going infeasible in 21% of
%                 hours against 5% at n = 4. Kept only for comparison.
%     'matched'   (default) blocks 2..m keep every unit's original local
%                 inflow EXCEPT the head, which absorbs -sum(mean local) so
%                 its net head inflow returns to block 1's 125.3 kcfs. Every
%                 block is then hydrologically identical to block 1: same
%                 head inflow, same within-block wetting profile, same tail
%                 release, full per-unit variance everywhere. This is the
%                 mode that makes the chain a fair replicate.
%
%   The intercept shift uses the stationary mean of the AR(1) local model,
%   mean = c/(1-phi), so c_new = c - delta*(1-phi) moves the mean by -delta
%   without touching the persistence or the innovation variance.
%
%   SHOCKS. Block b uses C.zseq circularly shifted by (b-1)*floor(T/m)
%   hours. A shift preserves each column's marginal exactly -- these are
%   conditioned, heavy-tailed sequences (WY2024 spans -9.75 to +15.19), not
%   Gaussian draws -- and decorrelates the blocks.
%
%   REQUIRES realtimeSimulation.m, a copy of realtimeSimulation.m with two changes.
%   Do not edit the original, which produced the headline runs.
%     line  29  n = 4  ->  n = numel(C.units)
%     line 109  the ARCH-X lookup keys archx_params.csv by unit name and
%               that file carries only the four real projects, so a replica
%               named MCN_2 returns an empty row. The copy strips a trailing
%               _<digits> before the lookup, so a replica inherits its base
%               unit's fitted parameters.

%% ------------------------------------------------------------ arguments
p = inputParser;
addParameter(p, 'LocalMode', 'zeromean', ...   % see note: 'matched' measured WORSE
    @(s) any(strcmpi(s,{'matched','zeromean','raw'})));
addParameter(p, 'Quiet', false);
parse(p, varargin{:});
localMode = lower(p.Results.LocalMode);
quiet     = p.Results.Quiet;

assert(isscalar(m) && m == round(m) && m >= 1, 'm must be a positive integer.');
n0 = numel(C.units);
T  = size(C.ref,1);
n  = n0*m;

perUnitT = {'ref','qLow','qHigh','obs','local','outObs','zseq'};
for f = perUnitT
    assert(isfield(C,f{1}), 'case is missing field %s', f{1});
    assert(size(C.(f{1}),2) == n0, '%s has %d columns, expected %d', ...
           f{1}, size(C.(f{1}),2), n0);
end
assert(size(C.ar,1)    == n0, 'ar must be %d x 3', n0);
assert(size(C.sigMid,1)== n0 && size(C.sigSd,1) == n0, 'sigMid/sigSd must be %d x 8', n0);
assert(numel(C.R)      == n0, 'R must have %d entries', n0);
assert(numel(C.sys)    == n0, 'sys must have %d entries', n0);

Ct = C;

%% -------------------------------------------------------- unit identity
Ct.units = cell(1,n);
for b = 1:m
    for i = 1:n0
        if b == 1
            Ct.units{(b-1)*n0+i} = C.units{i};
        else
            Ct.units{(b-1)*n0+i} = sprintf('%s_%d', C.units{i}, b);
        end
    end
end

% Positional, so every replica block is the same n0 units repeated. The
% suffixed names built above are now labels for row and variable names only;
% nothing looks a unit up by name any more.
Ct.sys = repmat(C.sys(:)', 1, m);

%% ------------------------------------------------------------- topology
Ct.upstream = [0, 1:(n-1)];
Ct.lag = repmat(C.lag(:)', 1, m);

%% ------------------------------------------------- per-unit time series
Ct.ref    = repmat(C.ref,    1, m);
Ct.qLow   = repmat(C.qLow,   1, m);
Ct.qHigh  = repmat(C.qHigh,  1, m);
Ct.obs    = repmat(C.obs,    1, m);
Ct.outObs = repmat(C.outObs, 1, m);
Ct.R      = repmat(C.R(:)',  1, m);

Ct.ar     = repmat(C.ar,     m, 1);
Ct.sigMid = repmat(C.sigMid, m, 1);
Ct.sigSd  = repmat(C.sigSd,  m, 1);

Ct.local  = repmat(C.local,  1, m);
meanLocal = mean(C.local, 1, 'omitnan');
blockGain = sum(meanLocal);                 % kcfs the block passes downstream
phiHead   = C.ar(1,2);                      % AR(1) persistence of the head unit

switch localMode
    case 'zeromean'
        for b = 2:m
            cols = (b-1)*n0 + (1:n0);
            Ct.local(:,cols) = C.local - meanLocal;
            Ct.ar(cols,1)    = 0;
        end
    case 'matched'
        for b = 2:m
            h = (b-1)*n0 + 1;               % head unit of this block
            Ct.local(:,h) = C.local(:,1) - blockGain;
            Ct.ar(h,1)    = C.ar(1,1) - blockGain*(1 - phiHead);
        end
    case 'raw'
        % nothing to do
end

%% ---------------------------------------------------------------- shocks
Ct.zseq = zeros(T, n);
shift = floor(T/m);
for b = 1:m
    cols = (b-1)*n0 + (1:n0);
    Ct.zseq(:,cols) = circshift(C.zseq, (b-1)*shift, 1);
end

if isfield(Ct,'Rcorr'), Ct = rmfield(Ct,'Rcorr'); end

Ct.tiledBlocks    = m;
Ct.tiledBaseUnits = C.units;
Ct.tiledLocalMode = localMode;
Ct.tiledBlockGain = blockGain;

%% ------------------------------------------------------------ diagnostic
if quiet, return; end

fprintf('tileCase: %d blocks x %d units = %d projects, LocalMode ''%s''\n', ...
        m, n0, n, localMode);
fprintf('  chain: %s\n', strjoin(Ct.units(1:min(n,8)), ' -> '));
if n > 8, fprintf('         ... -> %s\n', Ct.units{end}); end

headIn1  = mean(C.routedExo(1:T)) + meanLocal(1);
tailMean = mean(C.outObs(:,end), 'omitnan');
fprintf('  block 1: head inflow %.1f kcfs, block gain %.2f, tail %.1f kcfs\n', ...
        headIn1, blockGain, tailMean);
fprintf('  %-10s %14s %14s\n', 'block head', 'net head in', 'vs block 1');
for b = 1:m
    switch localMode
        case 'matched',  hin = tailMean - blockGain;
        case 'zeromean', hin = tailMean;
        case 'raw',      hin = tailMean + (b-2)*blockGain;
    end
    if b == 1, hin = headIn1; end
    fprintf('  %-10s %14.1f %+14.1f\n', Ct.units{(b-1)*n0+1}, hin, hin-headIn1);
end
if ~strcmp(localMode,'matched') && m > 1
    fprintf(['  NOTE LocalMode ''%s'' leaves downstream blocks wetter than block 1.\n' ...
             '       Check per-block storage against R before trusting any timing.\n'], ...
             localMode);
end

H = 24; NV = 10; NR = 10;
fprintf('  LP at H=%d: %d variables, %d rows (n=%d baseline: %d / %d)\n', ...
        H, n*H*NV, n*H*NR, n0, n0*H*NV, n0*H*NR);
end

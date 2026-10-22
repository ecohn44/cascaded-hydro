%% ========================================================================
% DATA LOADING AND CASE CONSTRUCTION
% ========================================================================

function C = dataload(root, year, lowerProbability, halfWindowDays)
% Load one water year and build its historical storage envelope
%
% Inputs:
%   root             - repository root; optional, see PATHS below
%   year             - water year to evaluate, e.g. 2024
%   lowerProbability - percentile for the lower envelope, 0.20 for P20
%   halfWindowDays   - half-width of the seasonal pooling window, in days
%
% Outputs:
%   C                - the case struct: unit parameters, topology, the
%                      hourly reference and envelope, observed storage and
%                      outflow, the fitted AR and sigma terms, and prices
%
% Called ONCE per sweep. The result is sliced across a parfor of hundreds of
% jobs, so nothing here may depend on cfg; that belongs in dispatchSetup.
%
% The reference and envelope are constructed separately for every evaluated
% water year. The evaluated year is left out, preventing its realized storage
% trajectory from leaking into its own operating target.
%
% At each hour of the water year:
%   1. Take the same hour of day within +/- 7 days of that seasonal date.
%   2. Pool those observations across the other seven water years.
%   3. Set the reference to their mean and the lower boundary to the
%      requested percentile (P20 for every published run).
%   4. Set the upper boundary to the physical usable-storage maximum.
%
% The input record contains WY2018-WY2025. For an evaluated water year, the
% other seven years form its retrospective leave-one-year-out training set.
%
% PATHS. root is optional. Resolution order:
%   1. the value passed in, if non-empty;
%   2. the environment variable DDU_DATA_ROOT, for a data tree held outside
%      the repository (a shared drive, or a scratch copy);
%   3. this file's own folder, which is the repository root.
% Rule 3 is the default so that a fresh clone runs with no configuration.

%% 1. Resolve the repository root
if nargin < 1 || isempty(root)
    root = getenv('DDU_DATA_ROOT');
    if isempty(root)
        root = fileparts(mfilename('fullpath'));
    end
end
assert(isfolder(fullfile(root,'data')), 'dataload:noData', ...
    ['No data/ folder under "%s". Pass the repository root explicitly, or ' ...
     'set DDU_DATA_ROOT.'], root);

%% 2. Read the input record
S = load(fullfile(root, 'data', 'inputs.mat'));

% conditioned_shocks.mat holds one standardised shock path per water year,
% which becomes C.zseq and is the realised forecast error the world delivers:
% realtimeSimulation evaluates e = WorldKappa * sigma * zseq(t,i). A sweep
% replaces it with a scenario from the paired bank (see sweepRunOne), so this
% default path is what a single direct call to realtimeSimulation uses.
Z = load(fullfile(root, 'data', 'conditioned_shocks.mat'));

%% 3. Case identity, topology and constants
C.year = year;
C.units = {'MCN','JDA','TDA','BON'};
C.upstream = [0 1 2 3];              % 0 means exogenous Priest Rapids
C.lag = [12 1 0 2];
C.KAF = 0.0826446;
C.CPOW = 0.08467 * 0.88;
C.spillPenalty = 0.5;
C.slackPenalty = 1e5;
C.zseq = Z.(sprintf('zseq_%d', year));

n = numel(C.units);
T = 8760;

%% 4. Per-unit parameters
% One struct-array element per unit, in C.units order, so a unit is addressed
% by position -- C.sys(i).umin -- exactly as every other per-unit quantity in
% C already is, and a whole column of them reads as one line:
%
%     umin = [C.sys.umin]';
%
% Fields, as stored in data/inputs.mat:
%   umin, umax     turbine release limits                    kcfs
%   rr_up, rr_dn   ramp limit per hour; rr_dn is NEGATIVE    kcfs/h
%   pmax           nameplate capacity                        MW
%   R              usable storage                            kaf
%   fb0, fb_slope  forebay elevation, affine in storage      ft, ft/kaf
%   tw_a, tw_b     tailwater elevation, affine in outflow    ft, ft/kcfs
%
% The input file also carries v0, W and sigma, which this controller never
% reads: initial storage comes from C.obs(1,:), and the forecast-error scale
% comes from C.sigMid / C.sigSd or from data/archx_params.csv. They are
% carried through untouched rather than dropped, so nothing downstream loses
% access to them.
%
% unitFields is a REQUIREMENT, not the full contents. Every name on it must be
% present, so a renamed field stops the run here with the unit named instead
% of surfacing as a missing-field error inside the per-hour LP build, while an
% added field is passed along without comment.
unitFields = {'umin','umax','rr_up','rr_dn','pmax','R', ...
              'fb0','fb_slope','tw_a','tw_b'};
rawParams = jsondecode(S.params);
fieldOrder = {};
for i = n:-1:1
    assert(isfield(rawParams, C.units{i}), 'dataload:noUnit', ...
        'data/inputs.mat carries no parameter block for unit "%s".', C.units{i});
    got = fieldnames(rawParams.(C.units{i}))';
    missing = setdiff(unitFields, got);
    assert(isempty(missing), 'dataload:unitFields', ...
        ['unit "%s" is missing parameter field(s) [%s]. dataload requires ' ...
         '[%s]; update unitFields here together with the call sites that ' ...
         'read them.'], C.units{i}, strjoin(missing, ', '), ...
        strjoin(unitFields, ', '));
    if isempty(fieldOrder)
        % Documented fields first, then whatever else the file carries, in
        % file order. Fixing one order for all units is what lets the struct
        % array be assigned element by element.
        fieldOrder = [unitFields, setdiff(got, unitFields, 'stable')];
    else
        assert(isempty(setxor(got, fieldOrder)), 'dataload:unitFieldsDiffer', ...
            ['unit "%s" carries a different field set from unit "%s". Every ' ...
             'unit block in data/inputs.mat must carry the same fields.'], ...
            C.units{i}, C.units{n});
    end
    C.sys(i) = orderfields(rawParams.(C.units{i}), fieldOrder);
end

%% 5. Per-unit hourly series, reference and envelope
C.ref = zeros(T,n); C.qLow = zeros(T,n); C.qHigh = zeros(T,n);
C.obs = zeros(T,n); C.local = zeros(T,n); C.outObs = zeros(T,n);
C.ar = zeros(n,3); C.sigMid = zeros(n,8); C.sigSd = zeros(n,8);
C.R = zeros(1,n);
availableYears = double(S.years(:)');
trainingYears = availableYears(availableYears ~= year);
C.trainingYears = trainingYears;
C.lowerEnvelopePercentile = 100*lowerProbability;
C.seasonalHalfWindowDays = halfWindowDays;
for i = 1:n
    u = C.units{i};
    project = C.sys(i);
    historicalStorage = zeros(T, numel(trainingYears));
    for j = 1:numel(trainingYears)
        historicalStorage(:,j) = col(S.(sprintf('%s_%d_obs',u,trainingYears(j))));
    end
    [C.ref(:,i), C.qLow(:,i)] = rolling_climatology( ...
        historicalStorage, lowerProbability, halfWindowDays);
    C.qHigh(:,i) = project.R * ones(T,1);
    C.obs(:,i) = col(S.(sprintf('%s_%d_obs',u,year)));
    C.local(:,i) = col(S.(sprintf('%s_%d_local',u,year)));
    C.outObs(:,i) = col(S.(sprintf('%s_%d_out_obs',u,year)));
    C.ar(i,:) = reshape(S.(sprintf('%s_%d_ar',u,year)),1,[]);
    C.sigMid(i,:) = reshape(S.(sprintf('%s_%d_sig_du_mid',u,year)),1,[]);
    C.sigSd(i,:) = reshape(S.(sprintf('%s_%d_sig_du_sd',u,year)),1,[]);
    C.R(i) = project.R;
end
%% 6. Clip the envelope, then the exogenous series and prices
C.qLow = max(C.qLow, 0);
C.qHigh = min(C.qHigh, C.R);
C.routedExo = col(S.(sprintf('MCN_%d_routed_exo',year)));
C.upReleaseExo = col(S.(sprintf('MCN_%d_up_release_exo',year)));
C.lmp = col(S.(sprintf('lmp_%d',year)));
end

function x = col(x)
x = double(x(:));
end

function q = row_quantile(X, probability)
% Linear-interpolated quantile along each row, toolbox-free
X = sort(X, 2);
n = size(X, 2);
position = 1 + (n-1)*probability;
lo = floor(position);
hi = ceil(position);
weight = position-lo;
q = (1-weight)*X(:,lo) + weight*X(:,hi);
end

function [reference, lower] = rolling_climatology(X, probability, halfWindowDays)
% Historical mean and quantile in a seasonal window
%
% X is hour x training-year. Reshaping to hour-of-day x day-of-water-year
% keeps diurnal patterns intact. The day window wraps around the start/end
% of the water year, so October 1 can use late-September observations.
nHours = size(X,1);
assert(mod(nHours,24)==0, 'Historical series must contain complete days.');
nDays = nHours/24;
storage = reshape(X, 24, nDays, size(X,2));
reference = zeros(nHours,1);
lower = zeros(nHours,1);

for day = 1:nDays
    windowDays = mod((day-halfWindowDays:day+halfWindowDays)-1,nDays)+1;
    for hour = 1:24
        samples = reshape(storage(hour,windowDays,:),1,[]);
        index = (day-1)*24+hour;
        reference(index) = mean(samples);
        lower(index) = row_quantile(samples, probability);
    end
end
end

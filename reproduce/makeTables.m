%% ========================================================================
% REPRODUCE: THE PAPER'S TABLES
% ========================================================================
%
% Aggregate the per-run summaries into one CSV per published table
%
%   Run it:  >> makeTables        (script, not a function)
%
% The summary files under results/runs/ hold ONE ROW PER RUN -- 96 rows for
% the headline tag. The paper's tables are aggregates of those rows:
% averaged over inflow scenarios, split by uncertainty framework, and
% converted to the units the paper reports. This script performs that
% reduction and writes the result to results/tables/, so every published
% table has one file that matches it column for column.
%
% Outputs, one per published table:
%
%   table_system_parameters.csv  Table II, the cascade as implemented
%   table_archx_fits.csv         Table III, conditional-variance fits
%   table_dispatch_epsilon.csv   Table IV, generation, violations, revenue
%                                and penalties against epsilon
%   table_benchmarks.csv         Table V, real-time dispatch benchmarks
%   table_scalability.csv        Table VI, solve cost against cascade size
%   PROVENANCE.md                for each file, its source and arithmetic
%
% Anything matching table_*.csv that a run does not write is deleted, so
% renaming or merging a table cannot leave a stale file behind.
%
% Requires no solver. Reads the distributed summary files, and calls dataload
% once for the parameter table, which takes a few seconds.

clear; clc;

ROOT = fileparts(fileparts(mfilename('fullpath')));
addpath(ROOT); addpath(fullfile(ROOT,'functions'));

PS   = fullfile(ROOT,'results','runs');
OUT  = fullfile(ROOT,'results','tables');
if ~isfolder(OUT), mkdir(OUT); end

MCM   = 1.23348;    % 10^6 m^3 per kaf (1 acre-foot = 1233.48 m^3)
THETA = 400;        % the published operating point
EPSB  = 0.05;       % the risk level the benchmark table is matched at
ZBON  = 2.2414;     % the Bonferroni multiplier at n = 4, eps = 0.05
HEAD  = 'BPA_HL2_TH0000400';

prov = {};          % provenance sections, written at the end
made = {};          % table files this run produced, used to prune stale ones

fprintf('writing the paper tables to results/tables/\n\n');

%% 1. Physical and forecast parameters
% Paper Table II. These are inputs rather than run output, so they are read
% from the case struct the controller itself builds: dataload with the
% published floor settings, which is what every reported run called. Reading
% them from data/inputs.mat by hand would let this file drift away from what
% was actually solved.
C = dataload(ROOT, 2024, 0.20, 7);
T0 = table(C.units(:), [C.sys.R]', [C.sys.pmax]', ...
           [C.sys.umin]', [C.sys.umax]', [C.sys.rr_dn]', [C.sys.rr_up]', ...
           C.lag(:), [C.sys.fb0]', [C.sys.fb_slope]', ...
           [C.sys.tw_a]', [C.sys.tw_b]', C.ar(:,1), C.ar(:,2), ...
    'VariableNames', {'unit','R_kaf','pmax_MW','umin_kcfs','umax_kcfs', ...
                      'rr_dn_kcfs_h','rr_up_kcfs_h','lag_h','fb0_ft', ...
                      'fb_slope_ft_kaf','tw_a_ft','tw_b_ft_kcfs', ...
                      'ar0_kcfs','ar1'});
writetable(T0, fullfile(OUT,'table_system_parameters.csv'));
made{end+1} = 'table_system_parameters.csv';
fprintf('  table_system_parameters.csv  %d rows\n', height(T0));
prov{end+1} = { ...
 '## table_system_parameters.csv'
 ''
 'Manuscript Table II. The four Lower Columbia projects as implemented,'
 'together with the head model and the mean inflow forecast.'
 ''
 '| Column | Meaning |'
 '|---|---|'
 '| `unit` | project: MCN, JDA, TDA, BON, in cascade order |'
 '| `R_kaf` | usable storage |'
 '| `pmax_MW` | nameplate capacity |'
 '| `umin_kcfs`, `umax_kcfs` | release limits, enforced as hard constraints |'
 '| `rr_dn_kcfs_h`, `rr_up_kcfs_h` | ramp limits, likewise hard |'
 '| `lag_h` | routing lag from the upstream project |'
 '| `fb0_ft`, `fb_slope_ft_kaf` | forebay elevation, affine in storage |'
 '| `tw_a_ft`, `tw_b_ft_kcfs` | tailwater elevation, affine in release |'
 '| `ar0_kcfs`, `ar1` | the AR(1) mean forecast on local inflow |'
 ''
 'Source: `dataload(ROOT, 2024, 0.20, 7)`, the same call every published run'
 'makes, reading `data/inputs.mat`. The routed upstream release enters the'
 'water balance explicitly rather than the forecast, so the mean model needs'
 'no separate release term. MCN heads the cascade, so its upstream inflow is'
 'the exogenous Priest Rapids record.'};

%% 2. Fitted conditional-variance coefficients
% Paper Table III. These are controller inputs rather than run outputs, so
% the file is the fit itself. Two columns are added: the paper reports omega
% and gamma SCALED BY THE UNCONDITIONAL VARIANCE, and reports one gamma that
% is the sum of the two exogenous coefficients. Without those columns the
% file cannot be read against the table.
af = fullfile(ROOT,'data','archx_params.csv');
if isfile(af)
    A = readtable(af);
    s2 = A.sigma_DIU.^2;
    A.omega_norm = A.omega ./ s2;
    A.gamma_norm = (A.gamma_release + A.gamma_ramp) ./ s2;
    writetable(A, fullfile(OUT,'table_archx_fits.csv'));
    made{end+1} = 'table_archx_fits.csv';
    fprintf('  table_archx_fits.csv         %d rows\n', height(A));
    prov{end+1} = { ...
     '## table_archx_fits.csv'
     ''
     'Manuscript Table III. The fitted conditional-variance coefficients,'
     'estimated by maximum likelihood on the AR(1) local-inflow residuals'
     'over eight water years.'
     ''
     '| Column | Meaning |'
     '|---|---|'
     '| `unit` | project: MCN, JDA, TDA, BON, in cascade order |'
     '| `omega`, `alpha` | constant and ARCH terms, in kcfs^2 and dimensionless |'
     '| `beta` | GARCH persistence; identically zero in every fit, so the variance is ARCH-X rather than GARCH-X |'
     '| `gamma_release`, `gamma_ramp` | coefficients on the upstream release and its ramp |'
     '| `sigma_DIU` | unconditional residual standard deviation, which is the constant the DIU framework uses |'
     '| `release_min`, `release_range`, `ramp_min`, `ramp_range` | min-max normalisation of the two exogenous predictors |'
     '| `n` | observations the fit was estimated on |'
     '| `omega_norm`, `gamma_norm` | the same fit in the manuscript''s scaling, added here |'
     ''
     'Source: `data/archx_params.csv`, read verbatim and extended with the'
     'two scaled columns below. This is a calibrated input distributed with'
     'the repository, estimated on the full eight-year record with no'
     'train/test split -- 70,080 hourly observations per project, which is'
     'the sample size the manuscript states. The estimation code is not part'
     'of this release, so nothing here refits it. `docs/fitting.md` gives the'
     'estimator and the normalisation, and relates these columns to the'
     'published table.'
     ''
     '### Reading this against Table III'
     ''
     'The manuscript reports the coefficients scaled by the unconditional'
     'variance, so that omega is read as a fraction of it, and reports a'
     'single exogenous coefficient rather than two:'
     ''
     '```'
     'omega_paper = omega / sigma_DIU^2'
     'gamma_paper = (gamma_release + gamma_ramp) / sigma_DIU^2'
     'alpha_paper = alpha                      (already dimensionless)'
     '```'
     ''
     'Those are the `omega_norm` and `gamma_norm` columns. The two exogenous'
     'predictors are min-max normalised onto the same scale before entering'
     'the variance, so summing their coefficients loses nothing. MCN heads'
     'the cascade and has no upstream decision, so both are zero there by'
     'construction.'};
end

%% 3. Generation, violations, revenue and penalties against epsilon
% Paper Table IV. Two rows per risk level, DIU above DDU, each averaged over
% the inflow scenarios. Earlier drafts split this across two tables; the
% columns are unchanged, so the figures below also reproduce those.
%
% Revenue is taken at the SCHEDULED tariff (revenue_sched_usd), not at the
% realised price (revenue_usd). The runs carry settlement noise, which shifts
% the realised figure by about 0.4 percent and is scoring-only; the published
% table reports the tariff. The penalty is the curtailed energy valued at the
% same tariff, using each run's own mean price.
f = fullfile(PS, HEAD, ['summary_' HEAD '.csv']);
assert(isfile(f), 'makeTables:noHeadline', 'missing %s', f);
H = readtable(f);
S = H(strcmpi(string(H.method),'SSH'), :);          % the table is SSH
eps = unique(S.eps);
pols = {'DIU','DDU'};

rows = {};
for q = 1:numel(eps)
    r = S(abs(S.eps - eps(q)) < 1e-12, :);
    for p = 1:numel(pols)
        a = r(strcmpi(string(r.policy), pols{p}), :);
        assert(~isempty(a), 'makeTables:noPolicy', ...
            'no %s rows at eps = %.3g in %s', pols{p}, eps(q), HEAD);
        curtail = mean(a.curtail_MWh);
        rows{end+1} = {eps(q), pols{p}, ...
                       mean(a.energy_MWh)/1e3, ...
                       mean(a.lower_ivi_kaf)*MCM, ...
                       mean(a.revenue_sched_usd)/1e6, ...
                       curtail*mean(a.price_mean)/1e3, ...
                       curtail};
    end
end
T1 = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'eps','model','gen_GWh','ivi_Mm3','revenue_Musd','penalty_kusd', ...
     'curtail_MWh'});
writetable(T1, fullfile(OUT,'table_dispatch_epsilon.csv'));
made{end+1} = 'table_dispatch_epsilon.csv';
fprintf('  table_dispatch_epsilon.csv   %d rows\n', height(T1));
prov{end+1} = { ...
 '## table_dispatch_epsilon.csv'
 ''
 'Manuscript Table IV. Generation, storage violation, revenue and'
 'curtailment penalty under the supporting-hyperplane algorithm at'
 'theta = 400. Two rows per risk level, one per uncertainty framework.'
 ''
 '| Column | Meaning |'
 '|---|---|'
 '| `eps` | chance-constraint risk level |'
 '| `model` | DIU or DDU forecast variance |'
 '| `gen_GWh` | annual generation |'
 '| `ivi_Mm3` | integrated violation below the P20 storage floor |'
 '| `revenue_Musd` | revenue at the scheduled BPA tariff |'
 '| `penalty_kusd` | value of the generation given up to recourse release cuts |'
 '| `curtail_MWh` | the curtailed energy the penalty is priced from |'
 ''
 ['Source: `results/runs/' HEAD '/summary_' HEAD '.csv`, SSH rows only,']
 'averaged over the four inflow scenarios.'
 ''
 '`curtail_MWh` is not a published column. It is written so that the'
 'penalty can be checked without re-reading the run summaries: the penalty'
 'is that energy priced at each run''s own mean tariff.'};

%% 4. Real-time dispatch benchmarks
% Paper Table V, at theta = 400 and eps = 0.05, every formulation matched to
% the same 1 - eps = 0.95 reliability.
%
% READ THIS BEFORE CHANGING THE POLICY FILTER. The benchmark rows carry the
% CONSTANT DIU forecast variance, because robust and scenario-based
% optimization fix their uncertainty sets and samples before the release
% schedule is known and so cannot express decision dependence at all. The
% chance-constrained formulation is the only one reported both ways, and the
% DDU row is the paper's result. Filtering these rows to DDU instead -- the
% obvious reading of a table captioned "benchmarks" -- shifts every column.
bf = fullfile(PS, 'benchmarks.csv');
if isfile(bf)
    B = readtable(bf);
    rules = {'none','robust','scenario'};
    names = {'Deterministic','Robust','Scenario-based'};
    rows  = {};
    for q = 1:numel(rules)
        r = B(strcmpi(string(B.rule), rules{q}) & ...
              strcmpi(string(B.policy),'DIU'), :);
        assert(~isempty(r), 'makeTables:noRule', ...
            'no DIU rows for rule "%s" in benchmarks.csv', rules{q});
        if strcmpi(rules{q},'none')
            pick = r;                               % z = 0, nothing to match
        else
            % Group by BOTH the set size and the discard count: the scenario
            % rule sweeps them together, and two groups can share a val while
            % reaching very different multipliers. The published row is the
            % group whose realised margin lands closest to the Bonferroni z,
            % which is how the paper holds reliability fixed across rules.
            [~, ~, gi] = unique([r.val, r.ndisc], 'rows');
            zByGroup   = accumarray(gi, r.margin_z, [], @mean);
            [~, b]     = min(abs(zByGroup - ZBON));
            pick = r(gi == b, :);
        end
        rows{end+1} = {names{q}, 'DIU', mean(pick.margin_z), ...
                       mean(pick.energy_MWh)/1e3, ...
                       mean(pick.lower_ivi_kaf)*MCM, ...
                       1e3*mean(pick.solve_mean_s)};
    end
    % The chance-constrained rows are the Bonferroni arm of the headline tag,
    % which the paper uses here so that all five rows share one scalar
    % multiplier rather than mixing SSH's cutting planes into the comparison.
    for p = 1:numel(pols)
        cc = H(strcmpi(string(H.method),'BON') & ...
               strcmpi(string(H.policy),pols{p}) & ...
               abs(H.eps - EPSB) < 1e-12, :);
        assert(~isempty(cc), 'makeTables:noBonferroni', ...
            'no BON + %s rows at eps = %.3g in %s', pols{p}, EPSB, HEAD);
        rows{end+1} = {'Chance-constrained', pols{p}, ZBON, ...
                       mean(cc.energy_MWh)/1e3, ...
                       mean(cc.lower_ivi_kaf)*MCM, ...
                       1e3*mean(cc.solve_mean_s)};
    end
    T2 = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'formulation','model','z','energy_GWh','ivi_Mm3','solve_ms'});
    writetable(T2, fullfile(OUT,'table_benchmarks.csv'));
    made{end+1} = 'table_benchmarks.csv';
    fprintf('  table_benchmarks.csv         %d rows\n', height(T2));
    prov{end+1} = { ...
     '## table_benchmarks.csv'
     ''
     'Manuscript Table V. Five ways of enforcing the storage bound, all at'
     'theta = 400 and eps = 0.05, matched to the same 1 - eps = 0.95'
     'reliability.'
     ''
     '| Column | Meaning |'
     '|---|---|'
     '| `formulation` | deterministic, robust, scenario-based or chance-constrained |'
     '| `model` | DIU or DDU forecast variance |'
     '| `z` | safety multiplier the formulation applies |'
     '| `energy_GWh` | annual generation |'
     '| `ivi_Mm3` | integrated violation below the floor |'
     '| `solve_ms` | mean solve time per 24-hour rolling horizon |'
     ''
     'Sources: `results/runs/benchmarks.csv` for the deterministic, robust'
     ['and scenario-based rows; the Bonferroni arm of `' HEAD '` for the']
     'two chance-constrained rows.'
     ''
     'The first four rows carry the constant DIU variance. Robust and'
     'scenario-based optimization fix their sets and samples before the'
     'release schedule is known, so they cannot express decision dependence;'
     'the chance-constrained formulation is the only one reported both ways,'
     'and its DDU row is the result the paper reports.'
     ''
     'For the robust and scenario rules the sweep covers several set sizes.'
     'The row written here is the one whose realised multiplier lands closest'
     'to the Bonferroni z = 2.2414, which is how the paper holds reliability'
     'fixed across formulations. `results/runs/benchmarks.csv` carries the'
     'whole sweep if you want to see the others.'};
else
    fprintf('  table_benchmarks.csv         SKIPPED (run reproduce/benchmarks)\n');
end

%% 5. Computational cost against cascade size
% Paper Table VI. Timing belongs to the machine that produced it, so nothing
% is distributed; the file appears once reproduce/scalability has been run.
sf = dir(fullfile(ROOT,'results','scalability*.csv'));
if ~isempty(sf)
    copyfile(fullfile(sf(1).folder, sf(1).name), ...
             fullfile(OUT,'table_scalability.csv'));
    made{end+1} = 'table_scalability.csv';
    fprintf('  table_scalability.csv        copied from %s\n', sf(1).name);
    prov{end+1} = { ...
     '## table_scalability.csv'
     ''
     'Manuscript Table VI. Solve cost against the number of units, from the'
     'synthetic 8- and 12-unit cascades that `functions/tileCase.m` builds by'
     'tiling the four calibrated projects.'
     ''
     ['Source: `results/' sf(1).name '`, produced by `reproduce/scalability`']
     'with FULLYEAR true, the default. The distributed copy is the run that'
     'the manuscript reports, carried over from the working tree rather than'
     'repeated here, because solve timings belong to the machine that'
     'produced them: an Apple M-series laptop, MATLAB R2025b, Gurobi 12,'
     'seven workers. Rerunning `reproduce/scalability` on other hardware'
     'replaces this file with that machine''s timings, which is the point of'
     'the table rather than a problem with it.'
     ''
     '| Column | Meaning |'
     '|---|---|'
     '| `n` | units in the cascade: 4 calibrated, 8 and 12 tiled |'
     '| `method` | BON or SSH |'
     '| `nvar`, `nrow` | decision variables and rows in the hourly program |'
     '| `attain_pct` | share of hours meeting the joint target exactly |'
     '| `t_mean`, `t_p95`, `t_max` | solve time per rolling horizon, seconds |'
     '| `iter_mean`, `iter_p95`, `iter_max` | SSH cutting-plane iterations; NaN for BON, which solves once |'
     '| `cuts_mean` | supporting hyperplanes added per horizon |'
     '| `wall_s` | wall clock for the whole year |'};
elseif isfile(fullfile(OUT,'table_scalability.csv'))
    % An earlier run on this machine produced it and the source CSV has since
    % been moved. Keep it rather than pruning a real timing result below.
    made{end+1} = 'table_scalability.csv';
    fprintf('  table_scalability.csv        kept from an earlier run\n');
else
    fprintf('  table_scalability.csv        SKIPPED (run reproduce/scalability)\n');
end

%% 6. Stale files
% Merging or renaming a published table leaves the old file behind, which is
% worse than no file at all: it still reads like a current result. Anything
% matching table_*.csv that this run did not write is removed.
old = dir(fullfile(OUT,'table_*.csv'));
for q = 1:numel(old)
    if ~ismember(old(q).name, made)
        delete(fullfile(OUT, old(q).name));
        fprintf('  removed stale %s\n', old(q).name);
    end
end

%% 7. Provenance
% Written last, so that it describes exactly the files this run produced.
% Every line goes through %s: the text carries per-cent signs and braces that
% fprintf would otherwise read as format specifiers.
fid = fopen(fullfile(OUT,'PROVENANCE.md'), 'w');
assert(fid > 0, 'makeTables:provenance', 'cannot write PROVENANCE.md');
head = { ...
 '# Paper tables'
 ''
 'One file per published table, aggregated from the per-run summaries in'
 '`results/runs/`. Regenerate at any time with:'
 ''
 '```matlab'
 'reproduce/makeTables        % seconds, no solver needed'
 '```'
 ''
 ['Generated by `reproduce/makeTables.m` on ' datestr(now,'yyyy-mm-dd') ...
  ', against manuscript TSTE-00904-2026.']
 ''
 'The source summaries hold one row per run -- 96 rows for the headline tag.'
 'These files are the averages over inflow scenarios that the paper reports.'
 ''
 '## Units'
 ''
 '| Quantity | Unit here | Converted from |'
 '|---|---|---|'
 '| Generation | GWh | `energy_MWh` / 1000 |'
 '| Storage violation (IVI) | 10^6 m^3 | `lower_ivi_kaf` x 1.23348 |'
 '| Revenue | M$ | `revenue_sched_usd` / 10^6 |'
 '| Curtailment penalty | k$ | `curtail_MWh` x `price_mean` / 1000 |'
 '| Solve time | ms | `solve_mean_s` x 1000 |'
 ''
 'The simulation works in thousand acre-feet. One acre-foot is 1233.48 m^3,'
 'so 1 kaf = 1.23348 x 10^6 m^3. Every volume in these files is 10^6 m^3,'
 'which is the unit the manuscript states, and no further conversion is'
 'needed to compare them.'
 ''
 'Every column is written at full double precision. The manuscript rounds to'
 'one decimal place, or to the nearest unit for the volumes, so round these'
 'figures before comparing them: 26566.9781594655 GWh is the published'
 '26,567.0, not a disagreement with it.'
 ''
 '## Which file is which table'
 ''
 '| Manuscript | File |'
 '|---|---|'
 '| Table II, the cascade as implemented | `table_system_parameters.csv` |'
 '| Table III, conditional-variance fits | `table_archx_fits.csv` |'
 '| Table IV, dispatch against the risk level | `table_dispatch_epsilon.csv` |'
 '| Table V, real-time dispatch benchmarks | `table_benchmarks.csv` |'
 '| Table VI, cost against cascade size | `table_scalability.csv` |'
 ''
 'Table I is the nomenclature and carries no data. The two parameter tables'
 'are inputs rather than run output; the other three are aggregates of'
 '`results/runs/`, which holds the raw sweeps every one of them reduces.'
 ''
 '---'};
fprintf(fid, '%s\n', head{:});
for q = 1:numel(prov)
    fprintf(fid, '\n');
    fprintf(fid, '%s\n', prov{q}{:});
end
tail = { ...
 ''
 '---'
 ''
 '## Calibrated inputs, not reconstructed here'
 ''
 'Tables II and III are inputs to every run rather than outputs of one. The'
 'AR(1) mean model lives in `data/inputs.mat` and the conditional-variance'
 'coefficients in `data/archx_params.csv`; both were estimated on the full'
 'eight-year record before any dispatch was solved, and both are read by the'
 'controller as distributed. The price signal is `data/bpa_pf24_rates.csv`,'
 'the published PF-24 rate applied without rescaling. None of the three is'
 'refitted by this script, and the estimation code is not part of this'
 'release.'};
fprintf(fid, '%s\n', tail{:});
fclose(fid);
fprintf('\n  PROVENANCE.md                written\n');
fprintf('\ndone. %d files in results/tables/\n', ...
        numel(dir(fullfile(OUT,'*')))-2);

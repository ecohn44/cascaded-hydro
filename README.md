# Cascaded hydropower dispatch under decision-dependent uncertainty

MATLAB implementation of a rolling-horizon dispatch controller for a cascade
of four Lower Columbia hydroelectric projects — McNary, John Day, The Dalles
and Bonneville — operating under a joint chance constraint on reservoir
storage. The inflow forecast error variance at each project is a function of
the controller's own upstream release decisions, so the uncertainty is
decision-dependent.

Accompanies IEEE Transactions on Sustainable Energy manuscript
TSTE-00904-2026. The repository contains the controller, the input data, the
summary results behind every published table, and the scripts that regenerate
them.

## Contents

1. [Formulation](#1-formulation)
2. [Installation](#2-installation)
3. [Usage](#3-usage)
4. [Function reference](#4-function-reference)
5. [Directories](#5-directories)
6. [Reproducing published tables](#6-reproducing-published-tables)
7. [Licence](#7-licence)

---

## 1. Formulation

A 24-hour look-ahead linear program is constructed and solved at each hour of
the water year. Only the first hour is implemented; the horizon then advances
and the program is rebuilt from the measured state.

### 1.1 Objective

Minimised over storage `V`, turbine release `u`, spill `sp`, power `p`, the
tracking deviation `d` and the storage slacks:

```
min  − Σ λ(t)·p(i,k)  +  c_sp·Σ sp(i,k)  +  Σ (θ/R_i)·d(i,k)  +  M·Σ slack(i,k)
```

with `d(i,k) = |V(i,k) − V_ref(i,k)|`, spill coefficient `c_sp = 0.5` and
slack penalty `M = 1e5`. The tracking weight is normalised by usable storage
`R_i`, which renders a single `θ` comparable across projects of differing
capacity. The published runs use `θ = 400`.

### 1.2 Constraints

Mass balance couples storage to local inflow and to upstream releases offset
by a per-project travel lag of up to 12 hours. Asymmetric ramp limits bound
the hour-to-hour change in release; release and power are bounded by nameplate
limits. Power is bilinear in storage and release, since head rises with
storage, and the product is replaced by its McCormick envelope.

### 1.3 Joint chance constraint

Storage is required to remain within its operating band at every project and
every look-ahead hour simultaneously:

```
P( V_min(i,k) ≤ V(i,k) ≤ V_max(i)   for all i, all k )  ≥  1 − ε
```

`V_min` is the seasonal P20 reference and `V_max` the physical usable
maximum. The two bounds differ in the recourse available against them, and
the implementation enforces this constraint in an equivalent form on the
lower bound alone: spill is unbounded above, so the upper bound holds for
every realisation of the inflow error and consumes no risk budget. The
reduction is exact rather than a relaxation. `docs/realtime.md` §8 derives
it, and `reproduce/twoSidedCheck` measures what enforcing both bounds costs.

Two enforcement methods are implemented.

**Bonferroni** divides the risk budget uniformly and displaces each bound by
`z·σ(i,k)`, with `z = Φ⁻¹(1 − ε/n)` over the `n` lower-bound tails. For
`n = 4` and `ε = 0.05`, `z = 2.2414`. The accumulated standard deviation
scales as `√k` over `k` hours.

**Sequential supporting hyperplanes** allocate the budget adaptively. From a
point strictly interior to the joint region, the method locates the
constraint boundary by bisection, evaluates the gradient of the joint
probability there, appends a supporting cut, and re-solves. Because `Σ` is
held fixed within a rolling step, the Bonferroni optimum satisfies every valid
cut, so the hyperplane objective cannot exceed it. Containment and objective
dominance are verified at every hour and reported in the summary.

### 1.4 Conditional variance

Forecast error variance follows an ARCH-X specification:

```
σ²(i,m) = ω_i(κ) + α_i·e²(i,m−1) + κ·( γᵘ_i·ũ(j) + γʳ_i·|Δũ(j)| )
```

where `ũ(j)` and `|Δũ(j)|` are the release level and ramp of the upstream
project `j`, min–max normalised over the estimation range. Within the
look-ahead, `e²` is unobserved for `m ≥ 2` and is replaced by the preceding
conditional variance, following the standard multi-step recursion.

`κ` scales the decision-dependent terms. `κ = 0` yields the
decision-independent control arm; `κ = 1` applies the fitted relationship. The
realised inflow error is generated from the fitted decision-dependent law in
both cases, evaluated at the controller's implemented ramps, so the
decision-independent arm is exposed to the dependence without modelling it.

Upstream ramps falling inside the look-ahead are predicted from the preceding
rolling plan, shifted one hour. That plan is data at time `t` rather than a
decision variable, which preserves convexity of the chance set and validity of
the cuts.

The complete formulation, with the derivations, is given in the
manuscript. `docs/realtime.md` is the corresponding code reference: what
each file builds, the variable ordering, and the reduction above.

---

## 2. Installation

| Requirement | Version | Required for |
|---|---|---|
| MATLAB | R2016b or later | all components |
| Gurobi | 13.0, MATLAB interface | all components |
| Parallel Computing Toolbox | any | `runSweep` with `Workers > 1` |
| Statistics and Machine Learning Toolbox | any | `quantile` in certain figure scripts only |
| Python | 3.9 or later | five figure scripts only |

Academic Gurobi licences are available at no cost from
<https://www.gurobi.com/academia/>.

```matlab
addpath('/Library/gurobi1303/macos_universal2/matlab')   % adjust to the local install
```

```bash
pip install -r docs/requirements.txt     # figure scripts only
```

## 3. Usage

Add the repository and its subfolders to the MATLAB path, then confirm the
installation:

```matlab
cd <repository root>
addpath(pwd); addpath(fullfile(pwd,'functions')); addpath(fullfile(pwd,'tests'))

checkEnvironment     % verifies toolchain, licence and input data
verifyParity         % rebuilds a published run and compares 56 fields
```

A single run is executed directly:

```matlab
C   = dataload(pwd, 2024, 0.20, 7);
cfg = struct('Horizon',24, 'Epsilon',0.05, 'Kappa',1, 'Theta',400, ...
             'MaxHours',168, 'SolutionMethod','SSH', ...
             'SigmaModel','archx', 'PowerModel','mccormick', ...
             'PriceMode','bpa', 'ChanceSides','lower');
out = realtimeSimulation(C, 'DDU', cfg);
disp(out.summary)
```

---

## 4. Functions

Signatures, call frequency and responsibility. The mechanism inside the
dispatch tiers -- the variable and row ordering, the margin rules, the
recourse ladder -- is in `docs/realtime.md`; the forecast and
conditional-variance models are in `docs/fitting.md`. This section is the
interface between them.

The controller is organised in three tiers, distinguished by call frequency.
Work is placed in the lowest-frequency tier that can perform it: a sweep runs
several hundred jobs over an 8760-hour year, so a quantity computed in the
wrong tier is recomputed by a factor of the job count or the hour count.

### 4.1 Pipeline

```
data/   inputs.mat   archx_params.csv   conditioned_shocks.mat
        bpa_pf24_lmp_wy2024.csv   scenario_bank_WY2024_*.mat
          |
          v
reproduce/<driver>        set cfg for one study, call runSweep     per study
          |
          v
runSweep                  enumerate the job grid                   per sweep
  dataload                build C, the water year                  per sweep
  parfor job
    sweepRunOne           run and save one job
      realtimeSimulation  own the state and the hour loop          per job
        dispatchSetup     build D, the job's fixed setup           per job
        hour loop         realtimeDispatch, solve, recourse        per hour
                          (docs/realtime.md)
      -> results/runs/<tag>/jobs/<job>.mat
          |
          v
collectResults            jobs -> summary_<tag>.csv                per tag
          |
          v
reproduce/makeTables      summaries -> results/tables/*.csv        per release
```

Every arrow crosses a file boundary, and every stage writes its output before
the next reads it, so a sweep can be interrupted and resumed at any stage.

### 4.2 Root functions

| Function | Called | Responsibility |
|---|---|---|
| `C = dataload(root, year, lowerProbability, halfWindowDays)` | once per sweep | Loads one water year and constructs its storage reference and envelope. Returns the case structure `C`: per-project parameters as a struct array `C.sys`, cascade topology and travel lags, hourly observed storage and outflow, the fitted autoregressive and variance terms, and prices. The reference and floor are built by leave-one-year-out seasonal climatology -- observations at the same hour of day within ±7 days are pooled across the other seven water years, the reference is their mean and the floor the requested percentile, P20 for every published run -- so the evaluated year is excluded from its own target. `root` is optional and resolves in order: the supplied value, the environment variable `DDU_DATA_ROOT`, then the folder containing `dataload.m`. |
| `D = dispatchSetup(C, policy, cfg)` | once per job | Builds everything fixed for the duration of one run: horizon and unit scaling, the flat variable index map, the constraint rows and bounds whose coefficients do not depend on the hour, the conditional-variance evaluator, the inflow-error correlation and the safety multiplier `z`. This is the only tier that reads `cfg`, and invalid configurations are rejected here rather than inside the hour loop. It is separate from `dataload` because `C` is built once and distributed across a `parfor` of several hundred jobs, while this setup depends on `cfg` and is rebuilt for each. |
| `out = realtimeSimulation(C, policy, cfg)` | once per job | Owns the reservoir state, the hour loop, the solve and the hyperplane handoff, the recourse ladder and the summary. It is the only function in which the simulation advances. `policy` is `'DIU'` or `'DDU'`. Returns `out.summary`, one structure of scalar metrics, and `out.trajectory`, the hourly series behind them. |
| `[bonModel, B] = realtimeDispatch(t, S, ctx)` | once per hour | Constructs that hour's linear program: projects local inflow over the look-ahead, evaluates the conditional standard deviation, displaces the storage bounds by the required margin, and assembles the model. `B` carries every intermediate the solve phase needs. The function is pure -- it reads state and context and returns a model, performs no solve, keeps no counters and carries nothing between hours -- which is what makes the per-hour build testable in isolation. |
| `[sg, s2] = forecastError(i, uLev, duR, s2prev, AX, omgK, kap, polK, sfloor, sceil)` | once per project per look-ahead step | Evaluates the conditional-variance model for one step. Returns the conditional standard deviation and the clipped variance to be fed to the next step of the recursion. |
| `[sol, model, diagOut] = applySSHNative(model, bonModel, gpar, ctx)` | once per hour, when `SolutionMethod` is `'SSH'` | Solves the rolling program under the joint chance constraint by supporting cuts. Candidate interior points are tried in order: the Bonferroni optimum, its max-centred variant from `findSlater`, then the max-centred point of the hyperplane model. If no candidate attains `1 - ε` at every look-ahead hour the joint chance set is empty from that state, as is the Bonferroni set it contains; the target is not relaxed, the same penalised-slack recourse solution is returned, and the hour is reported as `feasible = false` and counted. |
| `T = runSweep('Tag', tag, 'Epsilon', [...], 'Kappa', [...], ...)` | once per sweep | Enumerates a grid of jobs over scenario, ε, policy, κ and method, executes them in parallel and writes each to `results/runs/<Tag>/jobs/`. Finished jobs are skipped, so a sweep is resumable and may be extended by raising `Scenarios` and re-running the same tag. Settings must not be altered under an existing tag, since the folder would then hold two configurations that cannot be told apart. Inflow paths are drawn from a single shock bank at a fixed seed, so the two policy arms of a scenario face an identical realisation: comparisons are paired and reruns are bit-identical. |
| `T = collectResults(tag)` | once per tag | Reduces a finished sweep folder to `summary_<tag>.csv`, one row per run, and `paired_<tag>.csv`, one row per comparison. Settings are read from the `sweep_config.mat` written beside the jobs. |

### 4.3 Library

`functions/` holds the code shared between tiers. MATLAB local functions are
not visible outside their own file, so anything two files need lives here.

| Function | Called by | Responsibility |
|---|---|---|
| `model = makeModel(Aeq, Aub, beq, bub, obj, lb, ub)` | `realtimeDispatch`, `realtimeSimulation` | Assembles a Gurobi model structure from its constraint blocks. Pure: no state and no solve. |
| `[phi, gV, gSd, strain] = sshProbability(x, ctx, k)` | `applySSHNative`, `realtimeSimulation` | The joint storage reliability at look-ahead step `k` and its exact local sensitivities, by Genz separation of variables on a fixed Richtmyer lattice. Deterministic and smooth, which is what lets the hyperplane bisection converge. |
| `cache = sshProbabilityCache(Rcorr, numPoints)` | `dispatchSetup` | Precomputes the Cholesky factors and the quasi-Monte Carlo lattice that `sshProbability` needs, once per job rather than once per hour. |
| `[slater, mode] = findSlater(slaterModel, gpar, ctx)` | `applySSHNative` | A conservative point safely interior to the joint chance set, for the line search, used when the preselected minimum-release trajectory is too brittle to start from. |
| `[z, info] = scenarioMargin(n, Rcorr, S, r, seed)` | `dispatchSetup` | The scenario-program multiplier as an order statistic over `S` sampled paths with `r` discarded. Called with no arguments it prints the design table instead. |
| `Cs = sliceCase(C, h0, nHours)` | `reproduce/scalability`, `tests/` | Takes an `nHours` window of a case beginning at hour `h0`, leaving the unit set untouched. |
| `Ct = tileCase(C, m, ...)` | `reproduce/scalability`, `tests/` | Chains `m` copies of the calibrated case into a synthetic cascade of `4m` projects for the cost-against-size study. `LocalMode` selects how the replicas' local inflows are generated. |
| `sweepRunOne(job, C0, bank, cfg)` | `runSweep` | Runs one job and writes its result, substituting that job's shock path from the bank into the case. |

---

## 5. Directories

### `reproduce/`

One script per published result, each with its settings in a block at the
top, writing to its own tag and resumable. Listed with their outputs and run
times in §6. `analyzeBenchmarks` is the exception: it post-processes
`benchmarks` rather than producing a tag, comparing the benchmarks to the
hyperplane method at matched reliability.

### `tests/`

Five tiers, ordered by cost. `runAllTests` executes them in sequence and logs
to `logs/`.

| Script | Scope |
|---|---|
| `checkEnvironment` | Toolchain, Gurobi licence by trial program, input data, and every result file the manifest names |
| `testUnits` | Eleven invariants requiring no solver: the index map is a bijection, the Bonferroni quantile inverts the normal distribution, `κ = 0` removes decision dependence exactly, tiling replicates projects positionally, and every script resolves the repository root correctly |
| `verifyParity` | Rebuilds a published run and compares 56 summary fields. To be run after any change to the files in §4 |
| `smokeTestFullYear`, `smokeTestScale` | Behaviour and solve time on tiled cascades |
| `testSettlementNoise` | Confirms that settlement-price noise alters revenue without altering dispatch |
| `reproduceTables` | Re-runs a published table and compares it to the shipped summary, writing to `REPRO_<tag>/` |
| `paperTables` | The manifest: every result with its tag, driver, renderer and cost. Supplies `reproduceTables`; kept in step with section 6 |

### `functions/`

Code shared between tiers, tabulated with its callers in §4.3.

### `figures/`

Eight MATLAB and five Python scripts, writing PDF and PNG to `figures_out/`,
plus `variance_probe.py`, a module the two variance figures import for the
reconstruction and mean model so that neither refits anything. Most read the
summary files directly and require no solver. Three
(`make_tracking_figure`, `make_risk_allocation`, `window_figure_TDA`) require
hourly trajectory files that are not distributed; re-running the
corresponding driver regenerates them.

### `data/`

`inputs.mat` holds hourly inflow, outflow and forebay elevation for
WY2018–WY2025, obtained from USACE Dataquery, with the primary exports
retained in `data/raw/`, together with the physical parameters and the fitted
AR(1) mean model. `archx_params.csv` holds the ARCH-X conditional-variance
coefficients of Table III. Both are **calibrated inputs distributed with the
repository**: they were estimated on the full eight-year record, 70,080
hourly observations per project, and are read by the controller rather than
recomputed by it: no script in the release refits them. `docs/fitting.md`
gives the reconstruction, the two models, the estimator and the normalisation,
and relates the stored coefficients to the reported ones. Prices are the
Bonneville Power Administration PF-24
Load Shaping rate (Schedule PF-24 §2.1.3.1, FY2024–2025), applied as
published without rescaling.

### `fitting/`

`fitInflowFullYear.py` rebuilds the inflow reconstruction of
`docs/fitting.md` §1 and fits the conditional-variance specification of §3
under alternative samples and mean models, writing the candidate fits to
`fullyear_archx.csv`. Those rows are the specification comparison, not the
coefficients the controller reads, which are in `data/archx_params.csv`.

### `results/`

`tables/` holds one CSV per published table and is what to read; `runs/`
holds the raw sweeps those tables reduce, one folder per tag, each with
`summary_<tag>.csv` (one row per job) and the `sweep_config.mat` recording
the settings that produced it. `results/README.md` says what every tag is
for and which table or figure reads it.

### `docs/`

`realtime.md` is a code-level reference for the three dispatch tiers: what
each file builds, the variable and row ordering, the margin rules, and the
reduction of the joint chance constraint to its lower bound. The formulation
itself is not restated here; the manuscript carries it with the derivations.
`fitting.md` states the inflow mean model and the conditional-variance
specification as estimated, with the relation between the stored
coefficients and the published ones. `requirements.txt` pins the Python
versions the published figures were produced with; the MATLAB pipeline needs
none of them.

### `single-period/`

Archived work preceding this study: the single-period formulation, the earlier
Julia implementation, the streamflow model development, and the simulation
outputs behind them. It is retained for provenance and is not part of the
reproduction path -- nothing in the sections above reads from it, and
`single-period/README.md` describes what it holds and tabulates every absolute
path in it.

---

## 6. Reproducing published tables

### Table files

The summary files under `results/runs/` hold **one row per run**
-- 96 rows for the headline tag. The paper's tables are aggregates of those
rows. To produce one file per published table:

```matlab
reproduce/makeTables        % seconds, no solver
```

This writes to `results/tables/`:

| File | Manuscript table |
|---|---|
| `table_system_parameters.csv` | Table II, the four projects as implemented |
| `table_archx_fits.csv` | Table III, ARCH-X conditional-variance fits |
| `table_dispatch_epsilon.csv` | Table IV, generation, violations, revenue and penalties against ε |
| `table_benchmarks.csv` | Table V, real-time dispatch benchmarks at θ = 400, ε = 0.05 |
| `table_scalability.csv` | Table VI, computational cost against cascade size |
| `PROVENANCE.md` | For each file, its source rows, arithmetic and unit conversions |

Storage volumes convert as 1 kaf = 1.23348 x 10^6 m^3, so every volume in
these files is already the 10^6 m^3 the manuscript reports. Table I is the
nomenclature and carries no data. Tables II and III are calibrated inputs
rather than run output: both are built by reading what the controller itself
loads -- `dataload` with the published floor settings for the parameters, and
`data/archx_params.csv` for the variance coefficients -- so neither can drift
from what was actually solved, and neither is refitted here. The Table VI
timings are the distributed run from the authors' machine, since solve times
belong to the hardware that produced them.

### Source runs

The summary files are distributed with the repository; the scripts below
regenerate them from scratch. Run times are wall clock on seven workers.

| Result | File | Script | Time |
|---|---|---|---|
| Principal grid: DDU against DIU over six values of ε, both enforcement methods, θ = 400 | `results/runs/BPA_HL2_TH0000400/summary_BPA_HL2_TH0000400.csv` (96 rows) | `reproduce/headlineGrid` | 4.5 h |
| The same grid at θ = 100, 200 and 800 | `.../BPA_HL2_TH{0000100,0000200,0000800}/summary_*.csv` (96 rows each) | `reproduce/headlineGrid` | 13.5 h |
| Generation and violation against κ | `.../BPA_KAP_TH0000400_E050/summary_*.csv` (28 rows) | `reproduce/kappaSweep` | 1.6 h |
| Hours below target over θ = 200…1400, ten paired scenarios; the basis for θ = 400 | `.../BPA_TH7_P20_TH{0000200…0001400}/summary_*.csv` (20 rows each) | `reproduce/thetaSelection` | 1.8 h |
| Deterministic, robust and scenario model predictive control, one code path | `results/runs/benchmarks.csv` (100 rows) | `reproduce/benchmarks` | 1.5 h |
| Storage floor at P10 and P30 instead of P20 | `.../BPA_FLOOR_P{10,30}_TH0000400/summary_*.csv` (12 rows each) | `reproduce/floorSensitivity` | 20 min |
| Two-sided chance constraint, supporting §1.3 | `.../BPA_2SIDE_P10_TH{0000000…0000400}/summary_*.csv` | `reproduce/twoSidedCheck` | 10 min |
| Solve time and cut count at 4, 8 and 12 projects | `results/tables/table_scalability.csv` (6 rows) | `reproduce/scalability` | 2 h |

### Summary file columns

| Column | Definition |
|---|---|
| `scen`, `eps`, `policy`, `kappa`, `method` | Scenario index, ε, `DIU` or `DDU`, κ, `BON` or `SSH` |
| `energy_MWh`, `energy_adj_MWh` | Generation, and generation adjusted for terminal storage |
| `lower_ivi_kaf` | Integrated violation below the storage floor; the principal reliability metric |
| `lower_exit_pct`, `below_eps_hrs` | Proportion of hours below the floor; hours below the ε target |
| `rel_mean` | Mean realised joint reliability |
| `step_infeas` | Hours admitting no feasible plan from the current state |
| `revenue_usd`, `revenue_sched_usd` | Revenue at the realised and at the scheduled price |

### Automated comparison

```matlab
tests/reproduceTables('Only',{'BPA_HL2_TH0000400'})
```

The stored configuration for the named tag is replayed through the parser that
produced it, written to `REPRO_<tag>/`, and compared column by column against
the distributed summary. A correct rebuild is bit-identical: a single shock
bank at a fixed seed, one solver thread, a fixed solver seed and a linear
program at every step admit no drift. A rebuild that agrees only approximately
indicates that a non-deterministic element has been introduced.

---

## 7. Licence

MIT. See `LICENSE`.

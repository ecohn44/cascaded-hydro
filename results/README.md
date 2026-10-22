# Results

Two folders, with different jobs.

**`tables/`** is what to read. One CSV per published table, plus
`PROVENANCE.md` giving each file's source rows, arithmetic and units. Nothing
in here is raw: every file is already averaged over inflow scenarios and
converted to the units the manuscript prints. Regenerate the whole folder with
`reproduce/makeTables` (seconds, no solver).

**`runs/`** is what those tables are computed from, and it is needed. Each
folder is one sweep, named by the study that produced it, holding:

- `summary_<tag>.csv` — **one row per job**, a job being one (scenario,
  method, policy, risk level) combination run over the full water year, with
  39 recorded columns;
- `sweep_config.mat` — the exact settings that produced the tag, which
  `tests/reproduceTables.m` replays to check the shipped rows still come back.

Deleting `runs/` would leave the tables unverifiable and the figures
unbuildable: every figure script and `makeTables` itself reads these
summaries, and the reproduction tier has nothing to replay without the
configs. Together they are a few megabytes.

## Tag naming

`BPA_<study>_[<floor>]_TH<theta>` with theta zero-padded to seven digits, so
`TH0000400` is the published tracking weight of 400. `P10`, `P20` and `P30`
are the storage-floor percentile, and `E050` the risk level, where a study
varies them.

## What each tag is for

| Tag | Jobs | Study | Read by |
|---|---|---|---|
| `BPA_HL2_TH0000400` | 96 | the headline grid: DDU against DIU over the risk level, under both BON and SSH | **Tables IV and V**, `figures/make_ddu_heatmaps.m`, `figures/make_ssh_vs_bon.m` |
| `benchmarks.csv` | 28 | deterministic, robust and scenario-based MPC at matched reliability | **Table V**, `reproduce/analyzeBenchmarks.m` |
| `BPA_HL2_TH0000100`, `_TH0000200`, `_TH0000800` | 288 | tracking-weight stability either side of the published 400 | `figures/make_tracking_figure.m`, and theta 100 in `make_ddu_heatmaps.m` |
| `BPA_TH7_P20_TH0000200` … `_TH0001400` | 140 | tracking-weight selection across seven values of theta | `figures/make_fig_theta_select_hours.py` |
| `BPA_KAP_TH0000400_E050` | 28 | the decision-dependence sweep over kappa | `figures/make_kappa_figure.m` |
| `BPA_FLOOR_P10_TH0000400`, `_P30_` | 24 | floor-percentile sensitivity, P10 and P30 against the published P20 | read directly from the summaries |
| `BPA_2SIDE_P10_TH0000000` … `_TH0000400` | 8 | the two-sided chance constraint against the one-sided form | read directly from the summaries |

`tests/paperTables.m` is the machine-readable version of that table, and the
one `tests/checkEnvironment.m` checks against, so a tag that goes missing is
caught by tier 0 rather than by a figure failing months later.

## What is not here

Per-job hourly trajectories. One sweep can save them, and they run to
hundreds of megabytes, so they are not distributed;
`figures/window_figure_TDA.m` is the only script that reads them and it names
the missing file when run from a clone. Everything else — every table and
every other figure — builds from the summaries above.

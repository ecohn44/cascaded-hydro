# cascaded-hydro

MATLAB + YALMIP/Gurobi code for forward (non-anticipatory) simulation/optimization of a *cascaded* hydropower system under:
- **DET**: deterministic inflows
- **DIU**: decision-independent uncertainty (AR(1) mean + static covariance)
- **DDU**: decision-dependent uncertainty (AR mean + GARCH-X variance that depends on upstream releases)

The model supports two joint chance-constraint (JCC) approaches for storage (volume) reliability:
- **Bonferroni** (**jcc-bon**) static, pre-allocated risk budget
- **Supporting Hyperplane (SSH)** (**jcc-ssh**) adaptive risk attribution via cutting planes

---

## Optimization Repository layout

```
SOC-sims/
  simulation.m                # Main driver: data --> scenario --> optimization -->  plots
  dataload.m                  # System + seasonal forecasting parameters + scenario templates
  plotPolicies.m              # Plot/compare saved results (DET vs DIU vs DDU)
  plotEffectiveBounds.m       # Plot effective (chance-shifted) volume bounds
  resultsBonferroni/
    dry/  results_unit*_*.mat
    wet/  results_unit*_*.mat
  resultsSSH/
    dry/  results_unit*_*.mat
    wet/  results_unit*_*.mat

functions/
  genOptimization.m           # Core forward optimization loop (YALMIP + Gurobi)
  applySSH.m                  # SSH algorithm: reliability cuts + risk attribution weights
  findSlater.m                # Slater-point initializer for SSH feasibility
  scenarioSimulator.m         # Extreme-event inflow scenario generator (constant/pulse/extended)
  initSimSettings.m           # Validates and bundles sim settings
  simPlots.m, plotSSH.m, ...  # Plotting/reporting helpers
```

---

## Requirements

### MATLAB toolboxes
- **Optimization**: YALMIP (external) + **Gurobi** (solver)
- **Statistics and Machine Learning Toolbox**: for `mvncdf` and `norminv` (used in SSH / chance bounds)

### External dependencies
- **YALMIP** (https://yalmip.github.io/)
- **Gurobi** with MATLAB interface enabled

> In `simulation.m`, the first lines contain machine-specific `addpath(...)` calls.
> Update these to match your local install locations (YALMIP + Gurobi).

---

## Quick start

1. **Open MATLAB** and set your working directory to the `SOC-sims/` folder.
2. **Add paths** (either edit `simulation.m` or run these once per session):
   ```matlab
   addpath('<path-to-gurobi-matlab>');
   addpath(genpath('<path-to-yalmip>'));
   addpath(genpath(fullfile(pwd, '..', 'functions')));   % if functions/ sits next to SOC-sims/
   ```
3. **Run the main driver**:
   ```matlab
   simulation
   ```

You should see:
- Streamflow scenario plots
- A console 'VARIABLE REPORT'
- Trajectory plots for each unit (release/volume/etc.)
- If SSH is enabled: SSH reliability/alpha plots

---

## How to configure a run

All run configuration lives near the top of `simulation.m`.


### Core Simulation Settings
**The key settings to toggle between are** 
    % season - to switch between wet and dry seasns 
    % scenario - set as "pulse" for flood events and "extended" for drought events
    % uncertainty framework - to run uncertainty under different assumptions 
    % solution algorithm - switch between bonferroni "jcc-bon" and supporting hyperplane "jcc-ssh"

```matlab
simSettings = initSimSettings( ...
    "wet", ...        % season: "dry" | "wet"
    "pulse", ...      % scenario: "constant" | "pulse" | "extended" 
    "pwl", ...        % linear approximation method: "pwl" (head mapping via intervals)
    "diu", ...        % uncertainty framework: "det" | "diu" | "ddu"
    "jcc-bon", ...    % solution algorithm: "det" | "jcc-bon" | "jcc-ssh"
    "none" ...        % value of storage: "none" | "static" | "dynamic"
);
```


## What the model solves

Each time step `t = 1..T` solves a *single-step* LP (non-anticipatory) for all units:
- Decision variables per unit: storage `V`, power `p`, turbine release `u`, spill `s`
- Mass balance couples storage with inflow and outflows
- Ramp-rate limits enforce smooth changes in `u`
- Power uses a mapped hydraulic head `h(V)` (via interval mapping)
- Chance constraints shift storage bounds to enforce reliability (Bonferroni) or add SSH cuts

### Save toggles
```matlab
make_dir  = false;   % create timestamped plot folder
printplot = false;   % print/export plots (if implemented in plotting helpers)
save_mat  = false;   % save .mat results into resultsBonferroni/ or resultsSSH/
```
---

## Outputs

### In-memory outputs (from `genOptimization.m`)
- `X` (`T x 5n`): trajectories in blocks of 5 per unit
  - `V_i(t), p_i(t), u_i(t), s_i(t), sigma_i(t)`
- `std_hat` (`T x n`): predicted inflow standard deviation (DIU/DDU)
- `V_eff` (`T x 2n`): effective (chance-shifted) volume bounds
  - columns `[maxV_eff(i), minV_eff(i)]` per unit
- `phi_vals` (`T x 1`): achieved joint reliability under SSH (if `jcc-ssh`)
- `alpha_vals` (`T x n`): SSH risk attribution weights (sum to 1 across units)

### Saved results (`.mat`)
If `save_mat = true`, results are written to:
- `./resultsBonferroni/<season>/results_unit<k>_<framework>.mat`
- `./resultsSSH/<season>/results_unit<k>_<framework>.mat`

Saved variables include `X, V_eff, std_hat, q, sysparams, T, c, lag, season`.

---

## Plotting utilities

### Compare DET vs DIU vs DDU (pre-saved results)
Edit `season = "wet"` or `"dry"` in `plotPolicies.m`, and specify the solution alg to load files from then run:
```matlab
plotPolicies
```

### Effective bounds visualization
Specify the solution alg to load files from then run:
```matlab
plotEffectiveBounds
```


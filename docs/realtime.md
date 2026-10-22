# Real-time dispatch: code reference

Reference for the files that implement one rolling-horizon step. It covers
what each file builds, the variable and row ordering, the margin rules, and
the reduction of the joint chance constraint to its lower bound. The
formulation and its derivations are in the manuscript; this document is the
correspondence between that formulation and the code.

---

## Contents

1. [Call structure](#1-call-structure)
2. [`dataload`](#2-dataload)
3. [`dispatchSetup`](#3-dispatchsetup)
4. [`forecastError`](#4-forecasterror)
5. [`realtimeDispatch`](#5-realtimedispatch)
6. [`applySSHNative`](#6-applysshnative)
7. [`realtimeSimulation`](#7-realtimesimulation)
8. [Reduction to the lower bound](#8-reduction-to-the-lower-bound)
9. [Modifying this code](#9-modifying-this-code)

---

## 1. Call structure

Three tiers, separated by call frequency. A sweep runs hundreds of jobs over
an 8760-hour year, so work placed in the wrong tier is multiplied by the hour
count or the job count.

| Tier | File | Called | Returns |
|---|---|---|---|
| 1 | `dataload.m` | once per sweep | `C`, the water year |
| 2 | `dispatchSetup.m` | once per job | `D`, the job's fixed setup |
| 3 | `realtimeDispatch.m` | once per hour | the hour's linear program |

```
realtimeSimulation
    dispatchSetup                                   once, before the loop
    for t = 1:T
        realtimeDispatch(t, S, ctx)                 build
        gurobi(...) or applySSHNative(...)          solve
        recourse ladder, if required
        implement the first look-ahead step, advance the state
    assemble the summary
```

`realtimeDispatch` is pure: it reads state and context, returns a model, and
carries nothing between hours. All state resides in `realtimeSimulation`,
which is the only file that advances the simulation. The separation makes the
per-hour build testable in isolation.

---

## 2. `dataload`

```matlab
C = dataload(root, year, lowerProbability, halfWindowDays)
```

Loads one water year and constructs its storage envelope. Published runs call
`dataload(root, 2024, 0.20, 7)`.

`C` holds the unit parameters as a struct array `C.sys(i)`, the cascade
topology `C.upstream` and travel lags `C.lag`, the hourly storage reference
and envelope, observed storage and outflow, the fitted AR and variance terms,
and prices.

`dataload` is not folded into `dispatchSetup`. It is called once and the
resulting `C` is sliced across a `parfor` of hundreds of jobs, whereas the
setup depends on `cfg` and is rebuilt for each. Merging them would reload the
water year once per job.

---

## 3. `dispatchSetup`

```matlab
D = dispatchSetup(C, policy, cfg)
```

Builds everything that is fixed for the duration of one job: the variable
index map, the constraint rows whose coefficients do not depend on the
real-time hour, the safety multiplier, and the sigma evaluator.

`n` is taken from `numel(C.units)` rather than assumed. It is 4 for every
published run, but `functions/tileCase.m` constructs 8- and 12-unit cascades
for the scalability study, and a hardcoded 4 would solve only the first four
units.

### 3.1 Variable index map

One flat vector per rolling solve, addressed by a single map:

```matlab
ix = @(i,k,v) ((i-1)*H + (k-1))*NV + v;
nvar = n*H*NV;
```

`i` is the unit, `k` the look-ahead hour, `v` the variable code. Codes, in
order, with `NV = 10`:

| Code | Name | Quantity |
|---|---|---|
| 1 | `VOL` | end-of-hour storage, kaf |
| 2 | `REL` | turbine release, kaf/h |
| 3 | `SPL` | spill, kaf/h |
| 4 | `POW` | power, MW |
| 5 | `TRK` | absolute deviation from the storage reference |
| 6 | `RMP` | absolute hour-to-hour change in release |
| 7 | `SLO` | lower storage slack, penalised |
| 8 | `SHI` | upper storage slack, penalised |
| 9 | `MINSLK` | minimum-release slack, fixed at zero outside recourse |
| 10 | `WDEF` | forecast-water correction, fixed at zero outside recourse |

Under `PowerModel = 'mccormick'` a bilinear product variable is appended as
code 11 (`PROD`) and `NV` becomes 11.

Look-ahead hour `k = 1` is the hour implemented. Hours 2 to `H` exist so that
the controller sees the consequences of the current decision; they are
discarded and rebuilt at the next hour.

Tier 1 asserts that the map is a bijection onto `1..nvar` under both power
models. An off-by-one would overlay two variables and corrupt the program
without raising an error.

### 3.2 Constraint rows

One mass-balance equality per unit and look-ahead hour, followed by `NR`
inequality rows:

| Code | Name | Enforces |
|---|---|---|
| 1 | `RUP` | ramp-up limit |
| 2 | `RDN` | ramp-down limit |
| 3–4 | `RPOS`, `RNEG` | `RMP` as the absolute ramp |
| 5 | `RLO` | storage above the lower bound, less `SLO` |
| 6 | `RHI` | storage below the upper bound, less `SHI` |
| 7–8 | `RTP`, `RTN` | `TRK` as the absolute tracking deviation |
| 9 | `RPW` | power below the head-linearised limit |
| 10 | `RMIN` | release above the minimum, less `MINSLK` |
| 11–14 | `MC1`–`MC4` | McCormick envelope, `'mccormick'` only |

Rows whose coefficients are independent of the real-time hour are assembled
in section 8 of `dispatchSetup` and reused every hour. The head terms, the
ramp-reachable bounds and the chance margin depend on the current state and
are rebuilt per hour in `realtimeDispatch`.

The mass balance routes upstream water with the per-unit travel lag in
`C.lag`, which reaches 12 hours at McNary. The upstream contribution at
look-ahead index `k` is taken from source index `k-1-lag`; where that index is
negative the water is already in transit and is read from the implemented
history instead.

### 3.3 Margin rule

The safety multiplier is resolved once per job in section 6:

| `MarginRule` | `z` | Formulation |
|---|---|---|
| `chance` | `Phi^-1(1 - eps/nTails)` | the paper's controller |
| `none` | 0 | deterministic MPC |
| `robust` | `cfg.Gamma` | fixed uncertainty budget |
| `scenario` | an order statistic over sampled paths | Calafiore–Campi |

`nTails` is `n` under `ChanceSides = 'lower'` and `2n` under `'both'`, giving
`z = 2.2414` at `n = 4` and `eps = 0.05`. All four rules run through one code
path, which makes the benchmark table a comparison of margins rather than of
four controllers.

`robust` and `scenario` decouple `z` from `eps` by design. Neither carries a
calibrated probability. Deriving `Gamma` or the sample size from `eps` would
collapse both onto Bonferroni and remove the comparison.

---

## 4. `forecastError`

The conditional standard deviation of the inflow error:

```
sigma^2 = omega_i(kappa) + alpha_i * e^2_{i,m-1}
          + kappa * ( gu_i * ut_j + gr_i * |dut_j| )
```

`ut_j` and `|dut_j|` are the release level and ramp of the upstream project
`j`, min–max normalised over the estimation range and clipped at zero. Within
the look-ahead `e^2` is unobserved for `m >= 2` and is replaced by the
preceding conditional variance, following the standard multi-step recursion.

Upstream ramps falling inside the look-ahead are predicted from the preceding
rolling plan, shifted one hour. That plan is data at time `t` rather than a
decision variable, which preserves convexity of the chance set and validity
of the supporting cuts. Decision dependence still propagates between steps
through the implemented releases.

Two scaling parameters are distinct:

| | Scales | Values |
|---|---|---|
| `Kappa` | the controller's model of decision dependence | 0 is the DIU arm; 1 is the fitted relationship |
| `WorldKappa` | the realised world | 1 is the fitted model; above 1 is a misspecification test |

The DIU arm is evaluated inside the DDU world: the realised error follows the
fitted decision-dependent law in both arms, evaluated at each controller's own
implemented ramps. The DIU arm is exposed to the dependence without modelling
it.

Tier 1 asserts that at `Kappa = 0` the release level and ramp change sigma by
exactly zero, and that at `Kappa = 1` they do not. Were the parameter to stop
acting as a control, every DDU-against-DIU comparison would be confounded
without raising an error.

---

## 5. `realtimeDispatch`

```matlab
[bonModel, B] = realtimeDispatch(t, S, ctx)
```

Assembles the hour's program from the static rows in `D` and the current
state in `S`.

### 5.1 Objective

```matlab
obj(ix(i,k,POW)) = -lmp(tk);          % revenue, negated to minimise
obj(ix(i,k,SPL)) =  C.spillPenalty;   % spill, 0.5
obj(ix(i,k,TRK)) =  cfg.Theta / p.R;  % tracking, normalised by pool size
obj(ix(i,k,SLO)) =  C.slackPenalty;   % storage slack, 1e5
obj(ix(i,k,SHI)) =  C.slackPenalty;
```

Dividing the tracking weight by `p.R` makes `Theta` comparable across
projects of different pool size: a deviation is counted as a fraction of the
pool in which it occurs rather than in absolute kaf. The slack penalty is five
orders of magnitude above the price, so slack is used only where the program
would otherwise be infeasible.

`Theta = 400` is the calibrated value. README section 6 names the run that
selects it. The selection resolves on hours below target rather than on
integrated violation, since the hour count measures events while the integral
is dominated by individual severe inflow draws.

### 5.2 Chance margin

Predicted storage carries an accumulated standard deviation, and the lower
bound is displaced by `z` times it:

```matlab
sdConst(i,k) = KAF*sqrt(k)*sigDIU(i) + <ramp-dependent terms>
marginKnown  = z * sdConst(i,k)
```

The `sqrt(k)` is the accumulation of independent hourly errors in storage over
`k` hours.

---

## 6. `applySSHNative`

```matlab
[sol, model, diagOut] = applySSHNative(model, bonModel, gpar, ctx)
```

Bonferroni divides the risk budget uniformly. Sequential supporting
hyperplanes allocate it adaptively: from a point strictly interior to the
joint chance region, the method locates the constraint boundary by bisection,
evaluates the gradient of the joint probability there, appends a supporting
cut, and re-solves. The joint probability and its gradient are evaluated in
`functions/sshProbability.m` by Genz separation of variables on a fixed
Richtmyer lattice, which keeps the evaluation deterministic and smooth enough
for the bisection to converge.

Because `Sigma` is fixed within a rolling step, the Bonferroni optimum
satisfies every valid cut, so the hyperplane objective cannot be worse. A
dominance guard checks this at every hour and returns the Bonferroni solution
if numerics indicate otherwise; the expected count is zero and the realised
count is reported in the summary.

Where no candidate interior point attains `1 - eps` at every look-ahead hour,
the joint chance set is empty from that state, as is the Bonferroni set it
contains. The method does not relax its target: it returns the same
penalised-slack recourse solution Bonferroni uses and reports
`feasible = false`. Those hours are counted in the summary.

---

## 7. `realtimeSimulation`

```matlab
out = realtimeSimulation(C, policy, cfg)
```

Owns the state, the rolling loop, the recourse ladder and the summary. It is
the only file in which the simulation advances.

### 7.1 Recourse

Two physical mechanisms, applied after the solve, in order:

1. **Spill.** Water above the physical maximum is spilled. This is the
   mechanism that makes the upper storage bound hold pathwise.
2. **Curtailment.** Where realised storage would fall below the floor,
   release is cut back, bounded below by the ramp-reachable minimum
   `max(umin, uprev - |rr_dn|)`.

`rr_dn` is stored negative, so every use takes `abs()`. Without it the bound
lands above the current release and curtailment cannot fire; tier 1 asserts
the sign, because with the sign wrong the curtailment path is inactive and
violations are over-reported.

Beyond those two, an emergency ladder unlocks `MINSLK` and `WDEF` and
re-solves, only once the original program is proven infeasible. It minimises
the water correction first, then the release shortfall at that correction,
and only then the original objective at both minima. Each use is counted in
the summary column `emergency`; a healthy run reports zero.

---

## 8. Reduction to the lower bound

The model constrains storage to a band; the code constrains the lower bound
only. The two are the same inequality.

### 8.1 Statement

For units `i = 1..n` and look-ahead steps `t = 1..H`, with `V_min(i,t)` the
seasonal P20 reference and `V_max(i)` the physical usable maximum:

```
(P2)   P( V_min(i,t) <= V(i,t) <= V_max(i),   all i, all t )  >=  1 - eps
(P1)   P( V_min(i,t) <= V(i,t),               all i, all t )  >=  1 - eps
```

The claim is that (P1) and (P2) are equivalent. `ChanceSides` selects which
is built.

### 8.2 The upper bound holds for every realisation

Spill enters the program as `s(i,t) >= 0` with no upper bound, so for any
inflow error the recourse

```
s(i,t)  <-  s(i,t) + ( Vtilde(i,t) - V_max(i) )^+
```

is admissible, and yields

```
V(i,t) = min( Vtilde(i,t), V_max(i) )  <=  V_max(i)      for every realisation.
```

### 8.3 Equivalence

Let

```
A = { V(i,t) >= V_min(i,t)   for all i, t }
B = { V(i,t) <= V_max(i)     for all i, t }
```

Section 8.2 gives `P(B) = 1`, hence `P(Bᶜ) = 0`. Then

```
A ∩ Bᶜ ⊆ Bᶜ                  ⟹  P(A ∩ Bᶜ) = 0
P(A) = P(A ∩ B) + P(A ∩ Bᶜ)   =  P(A ∩ B)
```

so `P(A ∩ B) >= 1 - eps` if and only if `P(A) >= 1 - eps`, which is
(P2) ⟺ (P1). The reduction is an identity, not a relaxation.

### 8.4 Why the lower bound has no counterpart

Raising `V(i,t)` requires reducing `u(i,t) + s(i,t)`, and both are bounded
below:

```
s(i,t) >= 0
u(i,t) >= max( umin(i), u(i,t-1) - |rr_dn(i)| )
```

so the admissible upward correction is bounded by

```
s(i,t) + u(i,t) - max( umin(i), u(i,t-1) - |rr_dn(i)| )  <  infinity
```

and a sufficiently adverse realisation exceeds it. Hence `P(A) < 1` in
general while `P(B) = 1`. The asymmetry is physical: water can be released,
not created.

### 8.5 Cost of building (P2)

With `z(m) = Phi^-1(1 - eps/m)` the Bonferroni multiplier over `m` tails:

| `m` | Tails | `z` at `eps = 0.05` |
|---|---|---|
| `n = 4` | lower bounds | 2.2414 |
| `2n = 8` | both bounds | 2.4977 |

`z(2n)/z(n) = 1.114`. Every lower margin is inflated by 11.4% against an
event of probability zero, so (P2) as built defines a strictly smaller
feasible set and a weakly worse optimum than (P1).

### 8.6 Conditions

1. `V_max` must be the physical maximum. An operational envelope — a P80
   quantile, say — is not attainable by any admissible `s`, so `P(B) < 1`,
   section 8.2 fails, and the upper bound does require risk budget.
2. The reduction does not price spill at zero. Spilled water generates
   nothing and the objective charges `C.spillPenalty` for it, so the upper
   bound is carried as an economic cost rather than as a risk.

### 8.7 In the code

`ChanceSides = 'lower'`, the default and the published setting, sets
`nTails = n`. `'both'` sets `nTails = 2n` and builds (P2) directly, which is
what `reproduce/twoSidedCheck.m` runs. The summary records `spill_kaf` and
carries no upper-bound violation column, since section 8.2 makes that
quantity identically zero.

---

## 9. Modifying this code

Four invariants hold across the three tiers.

| | Invariant |
|---|---|
| Purity | `realtimeDispatch` reads state and returns a model. State belongs in `realtimeSimulation` and arrives through `S`. |
| Placement | A quantity fixed for the job is built in `dispatchSetup` and added to the pack in its section 9; a per-hour intermediate is added to `B` in section 3 of `realtimeDispatch` and unpacked in the caller. Over-returning is free, under-returning leaves a stale value in the caller. |
| Scope | Local functions are not visible outside their own file. A helper used by two files belongs in `functions/`. `checkcode` does not resolve cross-file calls and will not report a missing one. |
| Dimension | `n` is read from `numel(C.units)`. A literal silently truncates the tiled cascades built by `functions/tileCase.m`. |

Changes to these files are verified against the stored reference run, which
rebuilds one published job and compares 56 summary fields:

```matlab
clear functions
verifyParity
```

The tiers of the suite, what each establishes and what each costs are
documented in the header of `tests/runAllTests.m`.

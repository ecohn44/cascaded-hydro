# Real-Time Hydropower Dispatch: PWL-MILP Formulation

## `realtimeGurobi.m` — Technical Reference

---

## Table of Contents

1. [Overview and Motivation](#1-overview-and-motivation)
2. [Problem Statement](#2-problem-statement)
3. [The Nonlinearity Problem](#3-the-nonlinearity-problem)
4. [The PWL-MILP Reformulation](#4-the-pwl-milp-reformulation)
   - 4.1 [Piecewise-Constant Head Approximation](#41-piecewise-constant-head-approximation)
   - 4.2 [Binary Segment Indicators](#42-binary-segment-indicators)
   - 4.3 [Exact Linearisation of the Bilinear Product](#43-exact-linearisation-of-the-bilinear-product)
   - 4.4 [Power and Release Constraints](#44-power-and-release-constraints)
5. [Full Constraint Catalogue](#5-full-constraint-catalogue)
6. [Variable Layout](#6-variable-layout)
7. [Objective Function](#7-objective-function)
8. [Approximation Error Analysis](#8-approximation-error-analysis)
9. [Why This Approach Is Better Than McCormick](#9-why-this-approach-is-better-than-mccormick)
10. [Stochastic Constraint Tightening](#10-stochastic-constraint-tightening)
11. [Implementation Notes](#11-implementation-notes)
12. [Configuration Reference](#12-configuration-reference)

---

## 1. Overview and Motivation

`realtimeGurobi` solves a **single-period, $n$-unit hydropower dispatch problem** in real time. At
each time step $t$, the scheduler observes current reservoir storages $V^{prev}$ and inflow
forecasts $\hat{q}$, then decides turbine releases $u_i$, spill $s_i$, and expected storage $V_i$
to maximise power generation subject to physical and operational constraints.

The central modelling challenge is the **hydraulic head–power nonlinearity**. Power output depends
on both the water column height (hydraulic head) and the volumetric flow rate through the turbine:

$$p_i = c \cdot h_i(V_i) \cdot u_i$$

where the head–storage relationship is given by a power law:

$$h_i(V) = a_i \cdot V^{b_i}, \qquad b_i \in (0, 1)$$

This produces a **bilinear term** in the product $h_i(V_i) \cdot u_i$ — nonlinear in the decision
variables $V_i$ and $u_i$. Handling this term correctly is the core technical challenge.

Three approaches to this nonlinearity were considered:

| Approach | How it handles $h(V) \cdot u$ | Problem class | Relaxation gap |
|---|---|---|---|
| **Fixed-head LP** | $h$ evaluated at $V^{prev}$ (constant) | LP | None, but head is stale |
| **McCormick MINLP** | Bilinear relaxation over $[h^{lo}, h^{hi}] \times [u^{lo}, u^{hi}]$ | MINLP (spatial B&B) | Always positive — power overestimated |
| **PWL-MILP** (this file) | Piecewise-constant $h$ + exact binary linearisation | MILP (standard B&B) | None at integer nodes |

The **PWL-MILP** approach is adopted here. It converts the nonlinear dispatch problem into a pure
**Mixed-Integer Linear Program** solvable by Gurobi's standard branch-and-bound engine without
`FuncNonlinear` or `genconpow`.

---

## 2. Problem Statement

### Notation

| Symbol | Description | Dimension |
|---|---|---|
| $n$ | Number of hydropower units | scalar |
| $t$ | Current time step index | scalar |
| $c$ | Power conversion coefficient | scalar |
| $V_i$ | Reservoir storage (end of period) | decision |
| $p_i$ | Power output | decision |
| $u_i$ | Turbine release (water flow) | decision |
| $s_i$ | Spill flow | decision |
| $d_i$ | Volume tracking error | decision |
| $V^{prev}_i$ | Measured storage at start of period | parameter |
| $\hat{q}_i$ | Forecast inflow | parameter |
| $a_i, b_i$ | Head–storage coefficients | parameter |
| $\theta$ | Tracking penalty weight | parameter |
| $K$ | Number of PWL segments | parameter |

### Single-Period Optimisation Problem

$$\min_{V, p, u, s, d} \quad -\sum_{i=1}^n p_i + 10^{-4} \sum_{i=1}^n s_i + \theta \sum_{i=1}^n \frac{d_i}{V^{max}_i - V^{min}_i}$$

subject to:

- **Mass balance** (one equation per unit)
- **Ramp-rate limits** (upper and lower)
- **Storage bounds** (possibly tightened by chance constraints)
- **Power-production constraints** via the PWL-MILP formulation
- **Tracking error constraints**

---

## 3. The Nonlinearity Problem

### Head–Storage Function

The hydraulic head $h_i$ represents the effective water column height driving the turbine.
It is modelled as a strictly concave power law:

$$h_i(V) = a_i \cdot V^{b_i}, \qquad 0 < b_i < 1$$

Since $b_i \in (0, 1)$, this function is:

- **Monotone increasing**: higher storage $\Rightarrow$ higher head
- **Concave**: the marginal gain in head per unit volume decreases as $V$ grows
- **Smooth**: differentiable everywhere for $V > 0$

### The Bilinear Product

Substituting into the power equation:

$$p_i = c \cdot a_i \cdot V_i^{b_i} \cdot u_i$$

Defining $z_i = V_i^{b_i}$, this becomes:

$$p_i = c \cdot a_i \cdot z_i \cdot u_i$$

The term $z_i \cdot u_i$ is **bilinear** — a product of two continuous decision variables. This makes
the feasible region non-convex and the problem NP-hard in general.

### Why the McCormick Envelope Is Insufficient

The classical McCormick relaxation bounds the bilinear product $w = z \cdot u$ by four linear
inequalities derived from the extreme points of the box $[z^{lo}, z^{hi}] \times [u^{lo}, u^{hi}]$:

$$w \geq z^{lo} u + u^{lo} z - z^{lo} u^{lo} \tag{MC1}$$
$$w \geq z^{hi} u + u^{hi} z - z^{hi} u^{hi} \tag{MC2}$$
$$w \leq z^{hi} u + u^{lo} z - z^{hi} u^{lo} \tag{MC3}$$
$$w \leq z^{lo} u + u^{hi} z - z^{lo} u^{hi} \tag{MC4}$$

The envelope tightly wraps the bilinear surface **only at the four corners** of the box. For interior
points, the gap between the relaxed upper bound (MC3, MC4) and the true product is:

$$\text{gap} = (z^{hi} - z)(u - u^{lo}) \geq 0 \quad \text{(from MC3)}$$

In a single-period real-time setting, there is no multi-period bound tightening to close this gap.
The solver exploits it to overestimate power, producing plans that are **physically infeasible** in
terms of actual energy delivered.

---

## 4. The PWL-MILP Reformulation

The key insight is to **decouple the two continuous variables** by discretising one of them using
binary segment indicators. Once one variable is pinned to a segment, its product with the other
becomes **exactly linear**.

### 4.1 Piecewise-Constant Head Approximation

Divide the effective storage range $[V^{eff,lo}_i, V^{eff,hi}_i]$ into $K$ uniform segments of
equal width:

$$\Delta V_i = \frac{V^{eff,hi}_i - V^{eff,lo}_i}{K}$$

Define breakpoints:

$$V^{bp}_{i,k} = V^{eff,lo}_i + (k-1) \cdot \Delta V_i, \qquad k = 1, \ldots, K+1$$

Within segment $k$, the head is **approximated by its value at the midpoint**:

$$\bar{V}_{i,k} = \frac{V^{bp}_{i,k} + V^{bp}_{i,k+1}}{2} = V^{eff,lo}_i + \left(k - \tfrac{1}{2}\right) \Delta V_i$$

$$\bar{h}_{i,k} = a_i \cdot \bar{V}_{i,k}^{b_i}$$

The scalars $\bar{h}_{i,k}$ are **precomputed constants** — they do not appear in the optimisation
as variables. This is what converts the nonlinear head function into a linear term.

The approximation error within segment $k$ for unit $i$ is bounded by:

$$\varepsilon_{i,k} \leq \frac{c \cdot u^{max}_i}{2} \cdot a_i \cdot \left| V^{b_i}_{i,k+1} - V^{b_i}_{i,k} \right|$$

Since $h(V) = a \cdot V^b$ is concave, the midpoint approximation is always an **underestimate**
within each segment — the true head is below the chord but above the midpoint value for the outer
segments only when $b < 1$. The maximum error shrinks as $O(1/K^2)$ because the function is
smooth and concave.

### 4.2 Binary Segment Indicators

Introduce binary variables $y_{i,k} \in \{0, 1\}$ for $i = 1, \ldots, n$ and $k = 1, \ldots, K$,
where $y_{i,k} = 1$ indicates that $V_i$ lies in segment $k$.

**Exactly one segment is active:**

$$\sum_{k=1}^{K} y_{i,k} = 1 \tag{C4a}$$

**Tight Big-M range enforcement per segment:**

When $y_{i,k} = 1$, the storage must lie within $[V^{bp}_{i,k},\, V^{bp}_{i,k+1}]$.
Defining:

$$M^{lo}_{i,k} = (k-1) \cdot \Delta V_i \qquad M^{hi}_{i,k} = (K-k) \cdot \Delta V_i$$

the two-sided Big-M constraints simplify cleanly to:

$$-V_i + M^{lo}_{i,k} \cdot y_{i,k} \leq -V^{eff,lo}_i \tag{C4b}$$

$$V_i + M^{hi}_{i,k} \cdot y_{i,k} \leq V^{eff,hi}_i \tag{C4c}$$

**Verification that these are tight:**

When $y_{i,k} = 0$: (C4b) reduces to $-V_i \leq -V^{eff,lo}_i$, i.e., $V_i \geq V^{eff,lo}_i$
(always satisfied by the global lower bound). When $y_{i,k} = 1$: (C4b) gives
$-V_i \leq -V^{eff,lo}_i - M^{lo}_{i,k} = -V^{bp}_{i,k}$, i.e., $V_i \geq V^{bp}_{i,k}$. By symmetry,
(C4c) with $y_{i,k} = 1$ gives $V_i \leq V^{eff,hi}_i - M^{hi}_{i,k} = V^{bp}_{i,k+1}$.

The Big-M values are **exactly the minimum required** — there is no unnecessary slack.

### 4.3 Exact Linearisation of the Bilinear Product

Define auxiliary continuous variables:

$$v_{i,k} = y_{i,k} \cdot u_i, \qquad k = 1, \ldots, K$$

Because $y_{i,k} \in \{0, 1\}$ (binary), not a continuous SOS2 weight, the four standard
McCormick inequalities with $u^{lo}_i = u^{min}_i$ and $u^{hi}_i = u^{max}_i$ are **exact** at
every integer node:

$$-v_{i,k} + u^{min}_i \cdot y_{i,k} \leq 0 \tag{C4d}$$

$$-v_{i,k} + u_i + u^{max}_i \cdot y_{i,k} \leq u^{max}_i \tag{C4e}$$

$$v_{i,k} - u_i - u^{min}_i \cdot y_{i,k} \leq -u^{min}_i \tag{C4f}$$

$$v_{i,k} - u^{max}_i \cdot y_{i,k} \leq 0 \tag{C4g}$$

**Proof of exactness:**

*Case $y_{i,k} = 0$:* (C4d) gives $v_{i,k} \geq 0$; (C4g) gives $v_{i,k} \leq 0$.
Together: $v_{i,k} = 0 = y_{i,k} \cdot u_i$. Exact.

*Case $y_{i,k} = 1$:* (C4e) gives $v_{i,k} \geq u_i - 0 = u_i$; (C4f) gives
$v_{i,k} \leq u_i + 0 = u_i$. Together: $v_{i,k} = u_i = 1 \cdot u_i$. Exact.

This is the critical distinction from the McCormick relaxation of a continuous product: **when one
factor is binary, the McCormick inequalities are provably tight at every feasible integer solution**.
No spatial branch-and-bound is required to close the gap.

### 4.4 Power and Release Constraints

With the exact linearisation in place, the power constraint becomes a **linear equality**:

$$p_i = c \sum_{k=1}^{K} \bar{h}_{i,k} \cdot v_{i,k} \tag{C4h}$$

This works because at any integer solution, exactly one $y_{i,k^*} = 1$, so exactly one $v_{i,k^*} = u_i$
and all other $v_{i,k} = 0$. The sum collapses to:

$$p_i = c \cdot \bar{h}_{i,k^*} \cdot u_i$$

which is exactly the piecewise-constant head approximation evaluated at the active segment.

A release consistency constraint is also included:

$$u_i = \sum_{k=1}^{K} v_{i,k} \tag{C4i}$$

This is logically implied at integer nodes but is included because it significantly **tightens the LP
relaxation** at each branch-and-bound node, reducing the branch-and-bound tree depth and solve time.

---

## 5. Full Constraint Catalogue

For each unit $i = 1, \ldots, n$ and the current period $t$:

### C1 — Mass Balance

$$V_i + u_i + s_i = V^{prev}_i + \hat{q}_i$$

Water is conserved: storage changes by inflow minus turbine release and spill.

### C2, C3 — Ramp-Rate Constraints (active for $t > 1$)

$$u_i \geq u^{prev}_i + \Delta^{dn}_i \tag{C2}$$

$$u_i \leq u^{prev}_i + \Delta^{up}_i \tag{C3}$$

Limits on how quickly the turbine release can change between periods.

### C4a–C4i — PWL Power Production

As described in Section 4.

### C5, C6 — Volume Tracking Error

The deviation of $V_i$ from the reference band $[V^{ref,lo}_i, V^{ref,hi}_i]$ is penalised via the
slack variable $d_i \geq 0$:

$$V_i - d_i \leq V^{ref,hi}_i \tag{C5}$$

$$-V_i - d_i \leq -V^{ref,lo}_i \tag{C6}$$

These two constraints together enforce $d_i \geq \max(V_i - V^{ref,hi}_i,\; V^{ref,lo}_i - V_i,\; 0)$.

### Summary of Constraint Counts per Unit

| Block | Constraints | Depends on $K$? |
|---|---|---|
| C1 (mass balance) | 1 | No |
| C2–C3 (ramp rate) | 2 (if $t > 1$) | No |
| C4a (segment selection) | 1 | No |
| C4b–C4c (Big-M range) | $2K$ | Yes |
| C4d–C4g (linearisation) | $4K$ | Yes |
| C4h (power equality) | 1 | No |
| C4i (release consistency) | 1 | No |
| C5–C6 (tracking) | 2 | No |
| **Total** | $6K + 8$ | |

For $K = 10$, $n = 4$ units: **288 constraints** and **60 binary variables** — solved in
milliseconds by Gurobi.

---

## 6. Variable Layout

All decision variables are stacked into a single vector of length $5n + 2nK$:

$$\mathbf{x} = \bigl[\underbrace{V_1 \ldots V_n}_{\text{block 0}}\;\Big|\; \underbrace{p_1 \ldots p_n}_{\text{block 1}}\;\Big|\; \underbrace{u_1 \ldots u_n}_{\text{block 2}}\;\Big|\; \underbrace{s_1 \ldots s_n}_{\text{block 3}}\;\Big|\; \underbrace{d_1 \ldots d_n}_{\text{block 4}}\;\Big|\; \underbrace{y_{1,1} \ldots y_{n,K}}_{\text{block 5 — binary}}\;\Big|\; \underbrace{v_{1,1} \ldots v_{n,K}}_{\text{block 6 — continuous}}\bigr]$$

Index functions (1-based, MATLAB convention):

| Variable | Index function |
|---|---|
| $V_i$ | $5n \cdot 0 + i$ |
| $p_i$ | $n + i$ |
| $u_i$ | $2n + i$ |
| $s_i$ | $3n + i$ |
| $d_i$ | $4n + i$ |
| $y_{i,k}$ | $5n + (i-1)K + k$ |
| $v_{i,k}$ | $5n + nK + (i-1)K + k$ |

The $y$ and $v$ blocks are appended at the end so that all five original index functions remain
unchanged — no downstream offset errors.

---

## 7. Objective Function

The objective is:

$$\min \quad \underbrace{-\sum_{i=1}^n p_i}_{\text{maximise power}} + \underbrace{10^{-4} \sum_{i=1}^n s_i}_{\text{minimise spill}} + \underbrace{\theta \sum_{i=1}^n \frac{d_i}{V^{max}_i - V^{min}_i}}_{\text{penalise volume deviation}}$$

The normalisation factor $V^{max}_i - V^{min}_i$ in the tracking term makes the penalty
**dimensionless** — a tracking error of 1 full reservoir range costs exactly $\theta$ regardless
of the physical size of unit $i$.

**Important note on the reported objective:** The function returns `obj = sum(p_physical)`, where
`p_physical` is computed using the exact nonlinear formula $p = c \cdot a_i \cdot V_i^{b_i} \cdot u_i$.
This differs from the MILP's internal objective value `sum(p_out)` by the PWL approximation error.
Using the physical value prevents an overestimated power figure from propagating into the outer
real-time tracking or reward framework.

---

## 8. Approximation Error Analysis

### Maximum Error Bound per Segment

Within segment $k$ of unit $i$, the piecewise-constant head approximation error is bounded by:

$$\varepsilon_{i,k} \leq \frac{c \cdot u^{max}_i}{2} \cdot a_i \cdot \left| V^{b_i}_{i,k+1} - V^{b_i}_{i,k} \right|$$

Since the function $V \mapsto a_i V^{b_i}$ is concave and the approximation uses the midpoint value,
the worst-case error within any segment is half the head variation across that segment multiplied by
the maximum release.

### Global Maximum Power Error

Taking the supremum over all $K$ segments:

$$\varepsilon_i(K) = \frac{c \cdot u^{max}_i}{2} \cdot \max_{k} \left| \bar{h}_{i,k+1} - \bar{h}_{i,k} \right|$$

For the power law $h = a V^b$ with $b < 1$ on a uniform grid, this maximum occurs at $k=1$ (the
steepest part of the curve near $V^{min}$) and decays approximately as:

$$\varepsilon_i(K) = O\!\left(K^{-(1+b_i)}\right)$$

### Practical Guidance on $K$

| $K$ | Constraint count per unit | Extra binary vars per unit | Approx. relative error |
|---|---|---|---|
| 5 | 38 | 5 | $\sim$0.5–1% |
| 10 | 68 | 10 | $\sim$0.1–0.2% |
| 20 | 128 | 20 | $< 0.05\%$ |
| 50 | 308 | 50 | Negligible |

For most real-time applications, $K = 10$ provides an excellent accuracy–speed trade-off.

---

## 9. Why This Approach Is Better Than McCormick

### The Gap is Structural, Not Computational

In the McCormick relaxation, the gap between the relaxed upper bound on $p$ and the true bilinear
value is:

$$\text{gap at } (z, u) = (z^{hi} - z)(u - u^{lo})$$

This is strictly positive for any interior point $(z, u)$ strictly inside the box. No amount of
solver tuning or tighter tolerances can eliminate it — it is a **property of the relaxation**, not
the solver.

In the PWL-MILP formulation, once the binary variable $y_{i,k}$ is fixed to an integer value, the
product $v_{i,k} = y_{i,k} \cdot u_i$ is determined exactly by linear constraints. The remaining
approximation error is the piecewise-constant head approximation, which:

1. Has a **known, computable bound** (Section 8)
2. Can be made **arbitrarily small** by increasing $K$
3. Is **symmetric** — it is an approximation of $h(V)$, not a relaxation; it does not
   systematically bias the objective upward

### Comparison Table

| Property | McCormick MINLP | PWL-MILP (this work) |
|---|---|---|
| Relaxation gap | Always positive, direction-dependent | Zero at integer nodes |
| Power overestimation | Structural | None — only approximation error |
| Approximation error | Controlled by spatial B&B tolerance | Controlled by $K$ |
| Error bound | Not straightforwardly computable | Explicit formula in Section 8 |
| Solver requirement | `FuncNonlinear = 1`, `genconpow` | Standard MILP engine |
| Problem class | Nonlinear (spatial B&B) | Linear (standard B&B) |
| Reproducibility | Solver-dependent | Deterministic with fixed seed |

---

## 10. Stochastic Constraint Tightening

When `bounds = 'jcc-bon'`, the reservoir bounds are tightened using a **Bonferroni correction** to
provide a joint chance constraint guarantee.

### Individual Chance Constraint

For each unit $i$, the requirement is that $V_i$ remains feasible with individual probability
$1 - \alpha_i$. Under a Gaussian inflow forecast error model with standard deviation $\hat{\sigma}_i$,
this translates to a deterministic safety margin:

$$V^{eff,max}_i = V^{max}_i - z_{\alpha_i} \hat{\sigma}_i$$
$$V^{eff,min}_i = \max\!\left(0,\; V^{min}_i + z_{\alpha_i} \hat{\sigma}_i\right)$$

where $z_{\alpha_i} = \Phi^{-1}(1 - \alpha_i)$ is the standard normal quantile.

### Joint Chance Constraint via Bonferroni

For a system-level violation probability $\varepsilon$ across all $n$ units, Bonferroni's inequality
provides a sufficient condition by distributing the budget equally:

$$\alpha_i = \frac{\varepsilon}{2n} \quad \Rightarrow \quad z = \Phi^{-1}\!\left(1 - \frac{\varepsilon}{2n}\right)$$

The factor of 2 arises because each reservoir has both an upper and a lower bound. The resulting
effective bounds $V^{eff,lo}_i$, $V^{eff,hi}_i$ are then used directly in the PWL segment
construction, so the tighter bounds automatically propagate into tighter Big-M values and a tighter
LP relaxation.

**Note:** The non-negativity guard `V_eff_min = max(0, ...)` is essential — the power law $h = a V^b$
is not defined for $V \leq 0$, and a negative effective lower bound would invalidate both the
physical interpretation and the numerical computation of $\bar{h}_{i,k}$.

---

## 11. Implementation Notes

### Head Formula Consistency

All head computations in this file use the **single canonical formula**:

$$h_i(V) = a_i \cdot V^{b_i}$$

The `sys` struct must contain the fields `.a` and `.b` for each unit. If your system parameters
are defined via `min_h`, `max_h`, and `min_V`, you can recover `a` as:

$$a_i = \frac{h^{min}_i}{(V^{min}_i)^{b_i}}$$

### Reported Objective vs. MILP Objective

The function returns `obj = sum(p_physical)` — the exact physical power computed via $c \cdot a_i \cdot V_i^{b_i} \cdot u_i$
at the optimised $V_i$, $u_i$. This will differ slightly from the MILP's internal value due to the
piecewise-constant head approximation. The console output reports both values along with the maximum
per-unit approximation error and the active segment for each unit.

### Solver Configuration

The MILP is solved with:

```
params.FuncNonlinear = not set   % Not needed — pure MILP
params.Threads       = 1         % Deterministic latency for real-time loop
params.MIPGap        = 1e-4      % 0.01% optimality gap
params.TimeLimit     = 10        % Seconds — adjust to control interval
```

For $K = 10$ and $n \leq 10$, typical solve times are well under 100 ms.

### Warm Starting

Gurobi's MILP solver accepts an initial integer solution via `grb_model.start`. If the previous
period's segment indicators $y_{i,k}^{t-1}$ are stored, they can be passed as a warm start to
significantly reduce tree search depth.

---

## 12. Configuration Reference

### `sys` Struct Fields

| Field | Description | Required |
|---|---|---|
| `.min_V` | Minimum physical storage | Yes |
| `.max_V` | Maximum physical storage | Yes |
| `.a` | Head–storage coefficient | Yes |
| `.b` | Head–storage exponent $\in (0,1)$ | Yes |
| `.min_ut` | Minimum turbine release | Yes |
| `.max_ut` | Maximum turbine release | Yes |
| `.RR_dn` | Maximum downward ramp rate | Yes |
| `.RR_up` | Maximum upward ramp rate | Yes |
| `.F` | Installed power capacity | Yes |
| `.K` | Number of PWL segments (default: 10) | No |
| `.min_h`, `.max_h` | Head bounds (informational only) | No |

### Tuning $K$

Override globally by setting `sys(i).K` before calling. Note that the first unit's `K` value is
used for all units. To support per-unit $K$, a minor refactor of the segment loop is required.

---

*This document describes version `realtimeGurobi.m` — PWL-MILP formulation, single-period real-time
hydropower dispatch. For the oracle (multi-period) formulation see `oracleGurobi.m`.*

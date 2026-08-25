# `oracleNonlinear.m` — Technical Documentation

**Author:** Eliza Cohn  
**Version:** 1.0  
**Date:** August 2026  
**Solver:** Gurobi 13.0 MATLAB API  

---

## Table of Contents

1. [Purpose and Context](#1-purpose-and-context)
2. [Mathematical Formulation](#2-mathematical-formulation)
   - 2.1 [Decision Variables](#21-decision-variables)
   - 2.2 [Objective Function](#22-objective-function)
   - 2.3 [Constraints](#23-constraints)
   - 2.4 [Nonlinear Head-Storage Relationship](#24-nonlinear-head-storage-relationship)
3. [Function Signature](#3-function-signature)
   - 3.1 [Inputs](#31-inputs)
   - 3.2 [Outputs](#32-outputs)
4. [Variable Indexing](#4-variable-indexing)
5. [Constraint Inventory](#5-constraint-inventory)
   - 5.1 [Initial Conditions (t = 1)](#51-initial-conditions-t--1)
   - 5.2 [Dynamic Constraints (t >= 2)](#52-dynamic-constraints-t--2)
6. [Nonlinear Handling: genconpow and FuncNonlinear](#6-nonlinear-handling-genconpow-and-funcnonlinear)
7. [Gurobi Solver Parameters](#7-gurobi-solver-parameters)
8. [Solution Extraction and Output Matrix](#8-solution-extraction-and-output-matrix)
9. [Performance Benchmarks](#9-performance-benchmarks)
10. [Known Limitations and Assumptions](#10-known-limitations-and-assumptions)
11. [Usage Example](#11-usage-example)
12. [Architecture Notes: Role in the Broader Pipeline](#12-architecture-notes-role-in-the-broader-pipeline)

---

## 1. Purpose and Context

`oracleNonlinear` solves a **multi-unit hydropower dispatch optimization** over a finite
time horizon $T$. It computes the reference trajectory of reservoir volumes, power output,
water releases, and spill flows that maximize total energy generation subject to physical
and operational constraints.

The function replaces a YALMIP-based implementation that could not handle the nonlinear
head-storage relationship:

$$h_{i,t} = a_i \cdot V_{i,t}^{b_i}, \quad b_i \in (0, 1), \quad V_{i,t} > 0$$

YALMIP raises the error:

> `Warning: Solver not applicable (gurobi does not support signomial equality constraints)`

By calling the **Gurobi 13.0 MATLAB API directly**, the function uses `genconpow`
(power general constraints) combined with `FuncNonlinear=1` (spatial branch-and-bound)
to handle this nonlinearity exactly and globally.

This function serves as the **oracle** (reference trajectory generator) in a three-level
optimization pipeline:

```
oracleNonlinear (T-step reference)
    └── Real-time subproblem (per-step receding horizon)
         └── Supporting hyperplane algorithm (Joint Chance Constraints)
```

---

## 2. Mathematical Formulation

### 2.1 Decision Variables

| Symbol | Dimension | Description |
|--------|-----------|-------------|
| $V_{i,t}$ | $n \times T$ | Reservoir storage volume (per-unit) |
| $p_{i,t}$ | $n \times T$ | Power output at unit $i$, time $t$ |
| $u_{i,t}$ | $n \times T$ | Turbine water release (discharge) |
| $sp_{i,t}$ | $n \times T$ | Spill flow (excess release, penalized) |
| $z_{i,t}$ | $n \times (T-1)$ | Auxiliary variable: $z_{i,t} = V_{i,t}^{b_i}$, for $t \geq 2$ |

### 2.2 Objective Function

$$\max \sum_{i=1}^{n} \sum_{t=1}^{T} p_{i,t} - 10^{-4} \sum_{i=1}^{n} \sum_{t=1}^{T} sp_{i,t}$$

The spill penalty ($10^{-4}$) is small enough to not distort the primary generation
objective, but ensures spill is minimized as a tiebreaker between solutions of equal
energy value.

Internally Gurobi **minimizes**, so the model is assembled with objective coefficients
$-1$ on $p_{i,t}$ and $+10^{-4}$ on $sp_{i,t}$.

### 2.3 Constraints

**Bounds:**

$$s_i^{V,\min} \leq V_{i,t} \leq s_i^{V,\max} \quad \forall i, t$$

$$s_i^{u,\min} \leq u_{i,t} \leq s_i^{u,\max} \quad \forall i, t$$

$$0 \leq p_{i,t} \leq F_i \quad \forall i, t$$

$$sp_{i,t} \geq 0 \quad \forall i, t$$

**Initial conditions ($t = 1$):**

$$V_{i,1} + u_{i,1} + sp_{i,1} = V_{i,0} + I_{i,1}$$

$$u_{i,1} = u_i^{\min}$$

$$p_{i,1} = c \cdot h_{i,0} \cdot u_{i,1}, \quad h_{i,0} = a_i V_{i,0}^{b_i} \text{ (scalar constant)}$$

**Ramp rates ($t \geq 2$):**

$$\Delta^{\min}_i \leq u_{i,t} - u_{i,t-1} \leq \Delta^{\max}_i$$

**Mass balance ($t \geq 2$):**

For the headwater unit ($i = 1$) or when $t \leq \text{lag}$:

$$V_{i,t} - V_{i,t-1} + u_{i,t} + sp_{i,t} = I_{i,t}$$

For downstream units ($i \geq 2$, $t > \text{lag}$):

$$V_{i,t} - V_{i,t-1} + u_{i,t} + sp_{i,t} = I_{i,t} + u_{i-1,t-\text{lag}} + sp_{i-1,t-\text{lag}}$$

### 2.4 Nonlinear Head-Storage Relationship

The hydraulic head is modeled as:

$$h_{i,t} = a_i \cdot V_{i,t}^{b_i}$$

Since $b_i \in (0,1)$, this function is **strictly concave** and **increasing** in $V_{i,t}$.
Power production is:

$$p_{i,t} = c \cdot u_{i,t} \cdot h_{i,t} = c \cdot u_{i,t} \cdot a_i \cdot V_{i,t}^{b_i}$$

This is a bilinear product of $u_{i,t}$ and $h_{i,t}$, and nonlinear in $V_{i,t}$.

**Linearization strategy:** Introduce auxiliary variable $z_{i,t} = V_{i,t}^{b_i}$,
so $h_{i,t} = a_i z_{i,t}$. The power product becomes bilinear in $(u_{i,t}, z_{i,t})$,
handled by McCormick envelopes (see Section 5.2).

**McCormick Envelopes on $p_{i,t} = c \cdot a_i \cdot u_{i,t} \cdot z_{i,t}$:**

Let $\bar{h}^{\min}_i = a_i (V_i^{\min})^{b_i}$, $\bar{h}^{\max}_i = a_i (V_i^{\max})^{b_i}$,
$\bar{u}^{\min}_i$, $\bar{u}^{\max}_i$ be the known bounds. The four McCormick inequalities are:

$$\frac{p_{i,t}}{c} \geq \bar{h}^{\min}_i u_{i,t} + \bar{u}^{\min}_i a_i z_{i,t} - \bar{h}^{\min}_i \bar{u}^{\min}_i$$

$$\frac{p_{i,t}}{c} \geq \bar{h}^{\max}_i u_{i,t} + \bar{u}^{\max}_i a_i z_{i,t} - \bar{h}^{\max}_i \bar{u}^{\max}_i$$

$$\frac{p_{i,t}}{c} \leq \bar{h}^{\max}_i u_{i,t} + \bar{u}^{\min}_i a_i z_{i,t} - \bar{h}^{\max}_i \bar{u}^{\min}_i$$

$$\frac{p_{i,t}}{c} \leq \bar{h}^{\min}_i u_{i,t} + \bar{u}^{\max}_i a_i z_{i,t} - \bar{h}^{\min}_i \bar{u}^{\max}_i$$

These four constraints form the **convex hull** of the bilinear product over the box
$[\bar{h}^{\min}_i, \bar{h}^{\max}_i] \times [\bar{u}^{\min}_i, \bar{u}^{\max}_i]$,
providing the tightest possible linear relaxation of the power production term.

---

## 3. Function Signature

```matlab
[result, obj, X] = oracleNonlinear(T, c, I, lag, s)
```

### 3.1 Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `T` | `integer` | Number of time periods in the optimization horizon |
| `c` | `double` | Power conversion coefficient (scalar, unitless in per-unit space) |
| `I` | `double (n x T)` | Exogenous inflow matrix. `I(i,t)` is the natural inflow to unit $i$ at time $t$ |
| `lag` | `integer` | Cascade travel lag: number of time steps for water released from unit $i-1$ to reach unit $i$ |
| `s` | `struct array (1 x n)` | Per-unit system parameters (see table below) |

**Fields of `s(i)`:**

| Field | Description | Units |
|-------|-------------|-------|
| `s(i).a` | Head coefficient in $h = a V^b$ | per-unit |
| `s(i).b` | Head exponent, must satisfy $b \in (0,1)$ | dimensionless |
| `s(i).V0` | Initial reservoir storage volume | per-unit |
| `s(i).min_V` | Minimum allowable volume (must be $> 0$) | per-unit |
| `s(i).max_V` | Maximum allowable volume | per-unit |
| `s(i).min_ut` | Minimum turbine discharge | per-unit |
| `s(i).max_ut` | Maximum turbine discharge | per-unit |
| `s(i).RR_dn` | Ramp-down rate limit (negative value, e.g. $-0.1$) | per-unit per step |
| `s(i).RR_up` | Ramp-up rate limit (positive value, e.g. $+0.1$) | per-unit per step |
| `s(i).F` | Feeder (nameplate power) capacity | per-unit |

> **Critical requirement:** `s(i).min_V > 0` for all units. Gurobi's `genconpow` with
> a fractional exponent requires the argument variable to have a strictly positive lower
> bound. A value of `min_V = 0` will cause the solver to error or produce incorrect results.

### 3.2 Outputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `result` | `struct` | Full Gurobi result struct. Key fields: `result.status`, `result.x`, `result.objval` |
| `obj` | `double` | Total power generation over the horizon: $\sum_{i,t} p_{i,t}^*$ |
| `X` | `double (T x 5n)` | Solution matrix. Columns are ordered as $[V_1, p_1, u_1, sp_1, q_1, V_2, p_2, \ldots, q_n]$ where $q_i$ is the effective inflow (natural + cascade) |

**Column layout of `X`:**

| Columns | Content |
|---------|---------|
| $1$ to $5$ | $[V_1, p_1, u_1, sp_1, q_1]$ for unit 1 |
| $6$ to $10$ | $[V_2, p_2, u_2, sp_2, q_2]$ for unit 2 |
| $\vdots$ | $\vdots$ |
| $5(n-1)+1$ to $5n$ | $[V_n, p_n, u_n, sp_n, q_n]$ for unit $n$ |

If the solver does not find a feasible solution, `obj = NaN` and `X` remains empty.

---

## 4. Variable Indexing

All variables are stored in a single flat vector of length $N_{\text{vars}} = 4nT + n(T-1)$
using **0-based indexing** (required by the Gurobi MATLAB API).

```
Block 0 — V(i,t)  : indices  0            to  n*T - 1
Block 1 — p(i,t)  : indices  n*T          to  2*n*T - 1
Block 2 — u(i,t)  : indices  2*n*T        to  3*n*T - 1
Block 3 — sp(i,t) : indices  3*n*T        to  4*n*T - 1
Block 4 — z(i,t)  : indices  4*n*T        to  4*n*T + n*(T-1) - 1
          (defined only for t = 2..T)
```

**Index helper functions (anonymous, 0-based):**

```matlab
idx_V  = @(i,t)  (i-1)*T            + (t-1);
idx_p  = @(i,t)  n*T                + (i-1)*T + (t-1);
idx_u  = @(i,t)  2*n*T              + (i-1)*T + (t-1);
idx_sp = @(i,t)  3*n*T              + (i-1)*T + (t-1);
idx_z  = @(i,t)  4*n*T              + (i-1)*(T-1) + (t-2);  % t >= 2 only
```

> **MATLAB indexing note:** When addressing the `lb`, `ub`, `obj_coeff` arrays or the
> sparse matrix `A`, always add `+1` to convert from 0-based to MATLAB's 1-based indexing:
> e.g., `lb(idx_V(i,t)+1)`.

---

## 5. Constraint Inventory

### 5.1 Initial Conditions (t = 1)

Three constraints are added per unit $i$ at $t = 1$:

| # | Type | Expression | Form in Gurobi |
|---|------|-----------|----------------|
| 1 | Mass balance | $V_{i,1} + u_{i,1} + sp_{i,1} = V_{i,0} + I_{i,1}$ | `sense = '='` |
| 2 | Initial release | $u_{i,1} = u_i^{\min}$ | `sense = '='` |
| 3 | Power at $t=1$ | $p_{i,1} - c h_{i,0} u_{i,1} = 0$ | `sense = '='` |

Note: $h_{i,0} = a_i V_{i,0}^{b_i}$ is evaluated as a **scalar constant** (no decision
variable involved at $t=1$), so no `genconpow` constraint is needed at this time step.

### 5.2 Dynamic Constraints (t >= 2)

Per unit $i$ and time step $t \geq 2$, the following constraints are added:

| # | Type | Count per $(i,t)$ | Notes |
|---|------|-------------------|-------|
| 1 | Ramp-down bound | 1 | $u_{i,t} - u_{i,t-1} \geq \Delta_i^{\min}$ |
| 2 | Ramp-up bound | 1 | $u_{i,t} - u_{i,t-1} \leq \Delta_i^{\max}$ |
| 3 | McCormick LB 1 | 1 | $p/c \geq \bar{h}^{\min} u + \bar{u}^{\min} a z - \bar{h}^{\min}\bar{u}^{\min}$ |
| 4 | McCormick LB 2 | 1 | $p/c \geq \bar{h}^{\max} u + \bar{u}^{\max} a z - \bar{h}^{\max}\bar{u}^{\max}$ |
| 5 | McCormick UB 1 | 1 | $p/c \leq \bar{h}^{\max} u + \bar{u}^{\min} a z - \bar{h}^{\max}\bar{u}^{\min}$ |
| 6 | McCormick UB 2 | 1 | $p/c \leq \bar{h}^{\min} u + \bar{u}^{\max} a z - \bar{h}^{\min}\bar{u}^{\max}$ |
| 7 | Mass balance | 1 | See cascade logic below |
| **Total** | | **7 per $(i,t)$** | |

**Total linear constraint count:**

$$N_{\text{cons}} = 3n + 7n(T-1)$$

For $n = 6$, $T = 720$ (30 days): $N_{\text{cons}} = 18 + 30{,}114 = \mathbf{30{,}132}$

---

## 6. Nonlinear Handling: `genconpow` and `FuncNonlinear`

### `genconpow` Struct

For each unit $i$ and time step $t \geq 2$, a power general constraint is added:

```matlab
genconpow(k).xvar = idx_V(i,t);   % 0-based: argument variable V(i,t)
genconpow(k).yvar = idx_z(i,t);   % 0-based: result variable z(i,t)
genconpow(k).a    = s(i).b;       % exponent b_i in (0,1)
```

This instructs Gurobi to enforce $z_{i,t} = V_{i,t}^{b_i}$ exactly.

**Total `genconpow` constraints:** $n(T-1)$

For $n = 6$, $T = 720$: $6 \times 719 = \mathbf{4{,}314}$ power constraints.

### `FuncNonlinear = 1`

Setting this parameter activates Gurobi's **spatial branch-and-bound** algorithm:

- Instead of replacing $V^b$ with a static piecewise-linear grid before solving,
  Gurobi builds an **adaptive outer approximation** of $z = V^b$ dynamically at
  each node of the B&B tree.
- Cuts are added only where the relaxation is loose — near the optimal solution.
- This is both **exact** (no approximation error in the final solution) and
  **efficient** (few cuts needed in the unit domain $V \in [0,1]$).
- The final solution is **globally optimal** within the specified `MIPGap` tolerance.

> **Why this outperforms static PWL for this problem:**  
> In per-unit space $V \in [0,1]$ with $b = 0.45$, the curvature of $V^{0.45}$ is
> concentrated near $V = 0$. The spatial B&B naturally refines the approximation
> near the optimal $V^*$, whereas a static PWL grid wastes breakpoints in
> regions never visited by the solver.

---

## 7. Gurobi Solver Parameters

```matlab
params.OutputFlag    = 1;   % Print Gurobi log to console
params.Seed          = 1;   % Fix random seed for reproducibility
params.Threads       = 1;   % Single thread (set higher for production runs)
params.FuncNonlinear = 1;   % Spatial B&B for genconpow constraints
```

**Recommended additions for production use:**

```matlab
params.MIPGap        = 0.01;    % Relax from default 1e-4 to 1% (sufficient for oracle)
params.MIPFocus      = 1;       % Prioritise finding good primal solutions
params.TimeLimit     = 1800;    % 30-minute hard wall (seconds)
params.Threads       = 8;       % Use all available cores
params.NodefileStart = 0.5;     % Write B&B tree to disk if RAM > 0.5 GB
```

> **Note on `Threads = 1`:** The current setting enforces deterministic, reproducible
> solves. For production runs over long horizons (30–90 days), increasing to 4–8
> threads can reduce wall time significantly without affecting solution quality.

---

## 8. Solution Extraction and Output Matrix

After solving, the function checks `result.status` before extracting values:

```matlab
if strcmp(result.status, 'OPTIMAL') || strcmp(result.status, 'SUBOPTIMAL')
    % extract V_opt, p_opt, u_opt, sp_opt from result.x
else
    warning(...);
    obj = NaN;
end
```

**Effective inflow reconstruction:**

The output column `q(i,t)` represents the total effective inflow to unit $i$ at time $t$,
including natural inflow plus upstream cascade contributions:

$$q_{i,t} = \begin{cases} I_{i,t} & i = 1 \text{ or } t \leq \text{lag} \\ I_{i,t} + u_{i-1,t-\text{lag}}^* + sp_{i-1,t-\text{lag}}^* & i \geq 2,\; t > \text{lag} \end{cases}$$

This uses the **optimal release and spill values** from the upstream unit, so `q` reflects
the actual water availability at each downstream reservoir under the optimal policy.

---

## 9. Performance Benchmarks

The following benchmarks were obtained with `Threads = 1` on a standard workstation:

| Horizon | $n$ | $T$ | POW Constraints | B&B Nodes | Simplex Iters | Solve Time | MIPGap |
|---------|-----|-----|-----------------|-----------|---------------|------------|--------|
| 40 days | 6 | 960 | 5,754 | 2,356 | 102,985 | **9.71 s** | 0.0003% |

Key observations:
- The extremely tight gap (0.0003%) far exceeds the 1% requirement, suggesting
  `MIPGap = 0.01` could be used to terminate earlier and reduce solve time further.
- 2,356 nodes is a very small B&B tree for a problem of this size, confirming that
  the McCormick envelopes provide a tight root LP relaxation.
- Increasing `Threads` from 1 to 8 is expected to reduce wall time by 3–5x.

---

## 10. Known Limitations and Assumptions

| Limitation | Detail |
|------------|--------|
| `min_V > 0` required | Gurobi `genconpow` with fractional exponent requires strictly positive lower bound on the argument variable. Setting `min_V = 0` will cause incorrect behavior. |
| Single cascade chain | The mass balance assumes a single linear cascade ($i \to i+1$). Branching network topologies require structural changes to the mass balance constraints. |
| Deterministic inflows | Inflow `I` is treated as a known, deterministic matrix. Stochastic inflow handling (e.g., scenario trees) is not built into this function and is managed externally via the JCC pipeline. |
| Fixed lag | A single scalar `lag` applies uniformly to all cascade pairs $(i-1, i)$. Unit-specific lags would require per-pair indexing in the mass balance. |
| `b` exponent fixed | The head exponent `s(i).b` is a fixed scalar per unit. Time-varying exponents are not supported. |
| No integer variables | All variables are declared continuous (`vtype = 'C'`). Unit commitment (on/off) binary variables are not included. |

---

## 11. Usage Example

```matlab
%% Minimal test: 2 units, 10 time steps
n   = 2;
T   = 10;
c   = 0.9;
lag = 1;

for i = 1:n
    s(i).a     = 1.0;
    s(i).b     = 0.45;
    s(i).V0    = 0.5;
    s(i).min_V = 0.1;
    s(i).max_V = 1.0;
    s(i).min_ut = 0.05;
    s(i).max_ut = 0.50;
    s(i).RR_dn  = -0.10;
    s(i).RR_up  =  0.10;
    s(i).F      =  1.0;
end

I = 0.05 * ones(n, T);   % constant inflow

[result, obj, X] = oracleNonlinear(T, c, I, lag, s);

fprintf('Total generation: %.4f\n', obj);
fprintf('Solver status:    %s\n',   result.status);

%% Access per-unit results from X
% X has columns: [V1, p1, u1, sp1, q1, V2, p2, u2, sp2, q2, ...]
V1 = X(:, 1);   % Reservoir 1 volume trajectory
p1 = X(:, 2);   % Unit 1 power trajectory
```

---

## 12. Architecture Notes: Role in the Broader Pipeline

`oracleNonlinear` is the first of three optimization layers:

```
Layer 1 — Oracle (this function)
  Purpose : Generate reference trajectory (V*, u*, p*) over full horizon T
  Solver  : Gurobi spatial B&B, FuncNonlinear=1
  Output  : X matrix used to warm-start and inform layers 2 and 3

Layer 2 — Real-Time Subproblem (per time step)
  Purpose : Dispatch decision at each step t given current state
  Solver  : Same genconpow + FuncNonlinear=1, warm-started from Layer 1
  Note    : Only 6 POW constraints per solve — very fast

Layer 3 — Supporting Hyperplane Algorithm (JCC)
  Purpose : Enforce joint chance constraints on volume bounds
  Solver  : Pure LP master problem (no POW constraints)
  Note    : Gradient of p w.r.t. V is analytical:
            dp/dV = c * u * a * b * V^(b-1)
            evaluated using Layer 1 and 2 solutions as linearization points
```

**Gradient of power output with respect to reservoir volume** (for supporting hyperplane cuts):

$$\frac{\partial p_{i,t}}{\partial V_{i,t}} = c \cdot u_{i,t} \cdot a_i \cdot b_i \cdot V_{i,t}^{b_i - 1}$$

This is evaluated in closed form using the optimal $u_{i,t}^*$ and $V_{i,t}^*$ from the
oracle solution, requiring no additional nonlinear solves within the JCC iterations.

---

*End of documentation.*

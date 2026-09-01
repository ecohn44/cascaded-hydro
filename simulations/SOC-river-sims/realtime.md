# Real-Time Hydropower Dispatch Optimization

## Overview

This module implements a real-time, single-period hydropower dispatch
optimizer for a cascade of `n` hydro units over a horizon of `T`
periods. It uses a **McCormick envelope relaxation** of the nonlinear
power equation, solved as a pure **Linear Program (LP)** one period at
a time using Gurobi 13.0.

The system is organized as two cooperating layers:

```
PLANNING LAYER   oracleGurobi.m  (multi-period McCormick LP)
    |
    |  V_ref(i,t)  soft volume trajectory
    v
REAL-TIME LAYER  realtimeGurobi.m  (per-period McCormick LP)
    ^                   called by
    |
    realtimeSimulation.m  (simulation loop, state management)
```

Both layers use McCormick relaxations of the same nonlinear power
equation. The real-time layer achieves tighter relaxations because
`V_prev(i)` is a known measured constant at each call, allowing
the McCormick bounds to be tightened dynamically around the current
operating point.

---

## Files

| File | Role |
|---|---|
| `realtimeGurobi.m` | Single-period Gurobi LP oracle (McCormick) |
| `realtimeSimulation.m` | Outer simulation loop; calls oracle for t = 1,...,T |
| `oracleGurobi.m` | Multi-period McCormick LP planner (generates V_ref) |

---

## Mathematical Formulation

### Decision Variables (per period, per unit i)

| Variable | Symbol | Description |
|---|---|---|
| Storage | $V(i)$ | End-of-period reservoir volume |
| Power | $p(i)$ | Electrical power output |
| Release | $u(i)$ | Turbine water release |
| Spill | $\text{sp}(i)$ | Excess spill (non-negative) |
| Tracking error | $d(i)$ | Absolute deviation from reference volume |
| Head auxiliary | $z(i)$ | Linearised hydraulic head $\approx a_i \cdot V(i)^{b_i}$ |

### Objective

$$\min \quad -\sum_{i=1}^{n} p(i) \;+\; 10^{-4} \sum_{i=1}^{n} \text{sp}(i) \;+\; \theta \sum_{i=1}^{n} \frac{d(i)}{V^{\max}_i - V^{\min}_i}$$

The three terms represent:
1. Power maximisation (negative sign in minimisation objective)
2. Spill minimisation (light penalty to resolve degeneracy)
3. Volume tracking (soft guidance from the planning layer)

### Linear Constraints

| Label | Constraint | Type |
|---|---|---|
| C1: Mass balance | $V(i) + u(i) + \text{sp}(i) = V_{\text{prev}}(i) + I(i,t)$ | $=$ |
| C2: Ramp-down | $u(i) \geq u_{\text{prev}}(i) + \text{RR}_{\text{dn}}$ | $\geq$ |
| C3: Ramp-up | $u(i) \leq u_{\text{prev}}(i) + \text{RR}_{\text{up}}$ | $\leq$ |
| C4: Head linearisation | $z(i) = a_i b_i V_{\text{prev}}^{b_i-1} V(i) + a_i(1-b_i)V_{\text{prev}}^{b_i}$ | $=$ |
| MC1–MC4: McCormick | Four inequalities bounding $p(i) \approx c \cdot z(i) \cdot u(i)$ | $\leq$ |
| C5: Tracking upper | $V(i) - d(i) \leq V^{\text{ref}}(i)$ | $\leq$ |
| C6: Tracking lower | $-V(i) - d(i) \leq -V^{\text{ref}}(i)$ | $\leq$ |

Note that C1 is exact: $V_{\text{prev}}(i)$ is a known measured constant at
each call, not a decision variable. This is the central architectural
advantage of the real-time formulation.

### Power Equation and McCormick Relaxation

The true power output is:

$$p(i) = c \cdot a_i \cdot V(i)^{b_i} \cdot u(i)$$

This is nonlinear in two ways: the concave head term $V(i)^{b_i}$ and the
bilinear product with $u(i)$. Both are handled by introducing the
auxiliary variable $z(i) \approx a_i \cdot V(i)^{b_i}$ via a first-order
Taylor linearisation around $V_{\text{prev}}(i)$:

$$z(i) = a_i b_i V_{\text{prev}}^{b_i - 1} \cdot V(i) + a_i (1 - b_i) V_{\text{prev}}^{b_i}$$

The remaining bilinear product $p(i) = c \cdot z(i) \cdot u(i)$ is then
relaxed via four McCormick inequalities. Let:

$$z^L = a_i V^{L\,b_i}, \quad z^U = a_i V^{U\,b_i}, \quad u^L = \max(u^{\min}_i,\; u_{\text{prev}}(i) + \text{RR}_{\text{dn}}), \quad u^U = \min(u^{\max}_i,\; u_{\text{prev}}(i) + \text{RR}_{\text{up}})$$

where $V^L, V^U$ are tightened bounds around $V_{\text{prev}}$. The four
McCormick constraints are:

$$\frac{p(i)}{c} \geq z^L u(i) + u^L z(i) - z^L u^L $$

$$\frac{p(i)}{c} \geq z^U u(i) + u^U z(i) - z^U u^U $$

$$\frac{p(i)}{c} \leq z^U u(i) + u^L z(i) - z^U u^L $$

$$\frac{p(i)}{c} \leq z^L u(i) + u^U z(i) - z^L u^U $$

These four constraints form the tightest possible convex and concave
linear enclosure of the bilinear surface $z \cdot u$ over the box
$[z^L, z^U] \times [u^L, u^U]$.

### Bound Tightening Strategy

The key accuracy improvement over a static McCormick relaxation is
**dynamic bound tightening** at each period:

$$V^L_i = \max\!\left(V^{\min}_i,\; 0.8 \cdot V_{\text{prev}}(i)\right), \quad V^U_i = \min\!\left(V^{\max}_i,\; 1.2 \cdot V_{\text{prev}}(i)\right)$$

$$u^L_i = \max\!\left(u^{\min}_i,\; u_{\text{prev}}(i) + \text{RR}_{\text{dn}}\right), \quad u^U_i = \min\!\left(u^{\max}_i,\; u_{\text{prev}}(i) + \text{RR}_{\text{up}}\right)$$

Since $V_{\text{prev}}(i)$ is observed at each period, these bounds are
always centred on the true current state. The McCormick envelope width
shrinks proportionally, and the relaxation gap is correspondingly small.

---

## Convergence and Accuracy

The McCormick relaxation is exact when the bounds on $z$ and $u$ are
tight. The `PhysErr` diagnostic in `extractSolution` measures the
gap between the solver's $p(i)$ and the true physical power:

$$\text{PhysErr}(i) = p(i) - c \cdot a_i \cdot V(i)^{b_i} \cdot u(i)$$

A `PhysErr` of $\sim 0.1$ at system scale (4 units, $F_{\text{total}} \approx 15$)
corresponds to approximately a 0.7% relaxation gap, which is acceptable
for real-time dispatch.

### Comparison: Multi-Period vs. Real-Time McCormick

| Metric | Multi-period (90-day) | Real-time (per period) |
|---|---|---|
| Variables | $n \times T \times 6$ | $n \times 6$ |
| McCormick constraints | $4 \times n \times T$ | $4 \times n$ |
| Bound tightening | Static global $[V^{\min}, V^{\max}]$ | Dynamic around $V_{\text{prev}}$ |
| Model type | LP | LP |
| `NonConvex` required | No | No |

---

## Solver Parameters

| Parameter | Value | Reason |
|---|---|---|
| `MIPGap` | `1e-4` | 0.01% optimality gap for LP |
| `TimeLimit` | `10` | Hard per-period safety stop; should not trigger in normal operation |
| `Threads` | `0` | Use all available cores |
| `OutputFlag` | `0` | Suppress per-period Gurobi log |
| `Seed` | `1` | Fix random seed for reproducibility |

Note: `NonConvex = 2` and `FuncNonlinear = 1` are **no longer required**
because the model is a pure LP. These parameters have been removed from
`solveModel`.

---

## Initial Conditions

### V0 Initialisation

`V0` in `sysparams` must be set to match the historical reservoir levels
at the start of the simulation year, not a nominal fraction of capacity.
The long-term planner (`oracleGurobi.m`) generates `SOC_ref` from
historical data, and if `V0 << SOC_ref(:,1)`, the tracking penalty will
force large spills in the first period to try to close the gap.

The correct initialisation is to read `V0` from the first column of
`SOC_ref` (or the corresponding historical storage observation):

```matlab
for i = 1:n
    sysparams(i).V0 = SOC_ref(i, 1);
end
```

For the Columbia River cascade, the corrected values are approximately:

| Unit | Name | V0 (corrected) |
|---|---|---|
| 1 | McNary | 2.82 |
| 2 | John Day | 1.01 |
| 3 | Dalles | 0.86 |
| 4 | Bonneville | 1.46 |

---

## Usage

### Step 1: Generate the planning trajectory (McCormick planner)

```matlab
[result_plan, ~, X_plan] = oracleGurobi(T, c, I, SOC_init, theta, sys);

% Extract V_ref from the multi-period solution
V_ref = zeros(n, T);
for i = 1:n
    col_V = (i-1)*5 + 1;   % V is the first of 5 columns per unit in X
    V_ref(i,:) = X_plan(:, col_V)';
end
```

### Step 2: Align V0 with historical initial conditions

```matlab
for i = 1:n
    sysparams(i).V0 = SOC_ref(i, 1);
end
```

### Step 3: Run the real-time simulation

```matlab
results = realtimeSimulation(T, c, I, V_ref, SOC_init, theta, sysparams);
```

### Step 4: Access results

```matlab
results.p_history      % (n x T) power output
results.u_history      % (n x T) turbine release
results.V_history      % (n x T) end-of-period storage
results.sp_history     % (n x T) spill
results.p_physical     % (n x T) physically verified power
results.p_error        % (n x T) McCormick relaxation residual
results.track_error    % (n x T) normalised volume tracking error
results.runtime        % (1 x T) per-period Gurobi solve time (seconds)
results.failed_periods % list of t where fallback was used
```

---

## Solver Failure Fallback

If Gurobi returns a status other than `OPTIMAL` or `SUBOPTIMAL`,
`realtimeSimulation` applies a conservative fallback:

1. Hold the previous period release: $u(i,t) = u(i,t-1)$
2. Set spill to zero or to the minimum needed to respect $V^{\max}$
3. Compute $V(i,t)$ from the mass balance
4. Set $p(i,t) = 0$ (conservative; no power credited on fallback)
5. Log the period index in `results.failed_periods`

The state $V_{\text{prev}}$ still advances, so the next period's
oracle receives a valid measured state and can resume normal operation.

---

## Gurobi Version

Developed and tested on **Gurobi 13.0.3** (build v13.0.3rc0).

This formulation uses only standard LP constraints (`A`, `rhs`, `sense`,
`lb`, `ub`). No `genconstrnl`, `quadcon`, or nonlinear solver features
are used.

---

## References

- Gurobi Optimizer Reference Manual 13.0: Linear Programming
  https://docs.gurobi.com
- McCormick, G.P. (1976). Computability of global solutions to factorable
  nonconvex programs. *Mathematical Programming*, 10(1), 147–175.
- Gurobi Optimizer Reference Manual 13.0: MATLAB API
  https://docs.gurobi.com

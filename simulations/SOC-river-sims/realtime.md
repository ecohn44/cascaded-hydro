# Real-Time Hydropower Dispatch Optimization

## Overview

This module implements a real-time, single-period hydropower dispatch
optimizer for a cascade of `n` hydro units over a horizon of `T`
periods. It replaces the original multi-period McCormick relaxation with
an exact nonlinear formulation solved one period at a time using
Gurobi 13.0's spatial branch-and-bound engine.

The system is organized as two cooperating layers:

```
PLANNING LAYER   oracleGurobi.m  (multi-period McCormick LP)
    |
    |  V_ref(i,t)  soft volume trajectory
    v
REAL-TIME LAYER  realtimeGurobi.m  (per-period exact NLP)
    ^                   called by
    |
    realtimeSimulation.m  (simulation loop, state management)
```

---

## Files

| File | Role |
|---|---|
| `realtimeGurobi.m ` | Single-period Gurobi optimization oracle |
| `realtimeSimulation.m` | Outer simulation loop; calls oracle for t = 1,...,T |
| `oracleGurobi.m` | Multi-period McCormick planner (generates V_ref) |

---

## Mathematical Formulation

### Decision Variables (per period, per unit i)

| Variable | Symbol | Description |
|---|---|---|
| Storage | $V(i)$ | End-of-period reservoir volume |
| Power | $p(i)$ | Electrical power output |
| Release | $u(i)$ | Turbine water release |
| Spill | $\text{sp}(i)$ | Excess spill (non-negative) |
| Tracking error | $d(i)$ | Deviation from reference volume |
| Auxiliary | $w(i)$ | $V(i)^{b_i} \cdot u(i)$ (nonlinear auxiliary) |

### Objective

$$\min \quad -\sum_{i=1}^{n} p(i) \;+\; 10^{-4} \sum_{i=1}^{n} \text{sp}(i) \;+\; \theta \sum_{i=1}^{n} \frac{d(i)}{V^{\max}_i - V^{\min}_i}$$

The three terms represent:
1. Power maximization (negative sign in minimisation objective)
2. Spill minimization (light penalty to resolve degeneracy)
3. Volume tracking (soft guidance from the planning layer)

### Linear Constraints (6 per unit)

| Label | Constraint | Type |
|---|---|---|
| C1: Mass balance | $V(i) + u(i) + \text{sp}(i) = V_{\text{prev}}(i) + I(i,t)$ | $=$ |
| C2: Ramp-down | $u(i) \geq u_{\text{prev}}(i) - \text{RR}_{\text{dn}}$ | $\geq$ |
| C3: Ramp-up | $u(i) \leq u_{\text{prev}}(i) + \text{RR}_{\text{up}}$ | $\leq$ |
| C4: Power (linear in w) | $p(i) = c \cdot a_i \cdot w(i)$ | $=$ |
| C5: Tracking upper | $V(i) - d(i) \leq V^{\text{ref}}(i)$ | $\leq$ |
| C6: Tracking lower | $-V(i) - d(i) \leq -V^{\text{ref}}(i)$ | $\leq$ |

Note that C1 is exact: $V_{\text{prev}}(i)$ is a known measured constant at
each call, not a decision variable. This is the central architectural
advantage of the real-time formulation.

### Nonlinear Constraint (1 per unit, via `genconstrnl`)

$$w(i) = V(i)^{b_i} \cdot u(i) \qquad \forall\, i = 1, \ldots, n$$

This constraint is **non-convex** for two compounding reasons:

1. **Concavity of** $V^b$: for $b_i \in (0,1)$,

$$\frac{d^2}{dV^2} V^{b_i} = b_i(b_i - 1) V^{b_i - 2} < 0$$

2. **Indefinite Hessian of the bilinear product**: the Hessian of
$V \cdot u$ with respect to $(V, u)$ is:

$$H = [ 0 , 1 ; 1 , 0 ], \qquad \lambda_{1,2} = \pm 1$$

An indefinite Hessian is neither positive nor negative semi-definite,
confirming the constraint is neither convex nor concave. This is why
`NonConvex = 2` and `FuncNonlinear = 1` are both required in the
solver parameters.

---

## Expression Tree for `genconstrnl`

The constraint $w(i) = V(i)^{b_i} \cdot u(i)$ is encoded as a 5-node
directed tree in pre-order (parent before children):

```
         MULTIPLY          node 0  (root)
        /        \
      POW       VAR u(i)   nodes 1, 4
     /   \
VAR V(i) CONST b_i         nodes 2, 3
```

| Node | Opcode (int32) | Data (double) | Parent |
|---|---|---|---|
| 0 | 4 (MULTIPLY) | -1.0 | -1 (root) |
| 1 | 12 (POW) | -1.0 | 0 |
| 2 | 1 (VARIABLE) | `idx_V(i)` 0-based | 1 |
| 3 | 0 (CONSTANT) | $b_i$ | 1 |
| 4 | 1 (VARIABLE) | `idx_u(i)` 0-based | 0 |

**Critical indexing rule:**
- `resvar` uses a **1-based** MATLAB index
- `VARIABLE` node `data` uses a **0-based** Gurobi column index

---

## Convergence Analysis

Gurobi's spatial B&B constructs a linear outer
approximation (a McCormick-style polytope) of the nonlinear feasible
surface at each node. The tighter the bounds on $V(i)$ and $u(i)$, the
tighter this approximation is at the root node, and the fewer branching
steps are needed.

In the multi-period formulation, $V(i)$ ranges over its full global
domain $[V^{\min}_i, V^{\max}_i]$ at every time step, because future
storage is unknown. In the real-time formulation, $V_{\text{prev}}(i)$
is known, so the mass balance confines $V(i)$ to a narrow interval:

$$V^{\text{lb}}(i) = \max\bigl(V^{\min}_i,\; V_{\text{prev}}(i) + I(i,t) - u^{\max}_i\bigr)$$

$$V^{\text{ub}}(i) = \min\bigl(V^{\max}_i,\; V_{\text{prev}}(i) + I(i,t) - u^{\min}_i\bigr)$$

The width of this interval is at most $u^{\max}_i - u^{\min}_i$, which
is typically a small fraction of the total reservoir range. This means
the outer approximation of $V(i)^{b_i} \cdot u(i)$ is already very
tight at the root node, and the spatial B&B often terminates after
zero or very few branching steps.

### Comparison: Multi-Period vs. Real-Time

| Metric | Multi-period (90-day) | Real-time (per period) |
|---|---|---|
| Nonlinear constraints | $n \times (T-1) = 2876$ | $n$ |
| `quadcon` entries | $2876$ | $0$ |
| `NonConvex = 2` required | Yes | Yes |

---

## Solver Parameters

| Parameter | Value | Reason |
|---|---|---|
| `NonConvex` | `2` | Required to accept non-convex `genconstrnl` at model validation |
| `FuncNonlinear` | `1` | Activates spatial B&B for `genconstrnl`; always paired with `NonConvex = 2` |
| `MIPGap` | `1e-4` | 0.01% gap; achievable per period due to tight V bounds from V_prev |
| `TimeLimit` | `10` | Hard per-period safety stop; should not trigger in normal operation |
| `Threads` | `0` | Use all available cores; critical for per-period throughput |
| `OutputFlag` | `0` | Suppress per-period Gurobi log; simulation loop prints its own summary |
| `Seed` | `1` | Fix random seed for reproducibility of spatial B&B heuristics |

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

### Step 2: Run the real-time simulation

```matlab
results = runSimulation_realtime(T, c, I, V_ref, SOC_init, theta, sys);
```

### Step 3: Access results

```matlab
results.p_history     % (n x T) power output
results.u_history     % (n x T) turbine release
results.V_history     % (n x T) end-of-period storage
results.sp_history    % (n x T) spill
results.p_physical    % (n x T) physically verified power
results.p_error       % (n x T) solver relaxation residual (should be ~0)
results.track_error   % (n x T) normalised volume tracking error
results.runtime       % (1 x T) per-period Gurobi solve time (seconds)
results.mipgap        % (1 x T) per-period achieved MIPGap
results.failed_periods % list of t where fallback was used
```

---

## Solver Failure Fallback

If Gurobi returns a status other than `OPTIMAL`, `SUBOPTIMAL`, or
`TIME_LIMIT` with an incumbent, `runSimulation_realtime` applies a
conservative fallback:

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

`genconstrnl` is the preferred nonlinear constraint interface as of
Gurobi 13.0. The previously used `genconpow` interface is deprecated
in this version and should not be used in new code.

---

## References

- Gurobi Optimizer Reference Manual 13.0: Nonlinear Constraints
  https://docs.gurobi.com
- Gurobi Optimizer Reference Manual 13.0: Expression Trees
  https://docs.gurobi.com
- Gurobi Optimizer Reference Manual 13.0: genconstrnl (MATLAB API)
  https://docs.gurobi.com

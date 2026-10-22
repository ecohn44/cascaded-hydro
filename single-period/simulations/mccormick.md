# Real-Time Hydropower Dispatch: McCormick Envelope Formulation

## `realtimeGurobiMC.m`

---

## Overview

`realtimeGurobiMC` solves a **single-period, $n$-unit hydropower dispatch** at each time step $t$.
Given measured reservoir storages $V^{prev}$ and a forecast inflow $\hat{q}$, it decides turbine
releases $u_i$, spill $s_i$, and end-of-period storage $V_i$ to maximise physical power generation.

This file is a **testing and benchmarking variant**. It applies a McCormick envelope relaxation
to the nonlinear power production term. It is intentionally simple and is intended for direct
comparison against the PWL-MILP formulation in `realtimeGurobiPWL.m`.

---

## The Nonlinearity

Power output depends on hydraulic head and turbine release:

$$p_i = c \cdot h_i(V_i) \cdot u_i$$

where the head-storage relationship is:

$$h_i(V) = a_i \cdot V^{b_i}, \qquad b_i \in (0, 1)$$

Introducing $z_i = V_i^{b_i}$, the power term becomes a **bilinear product**:

$$p_i = c \cdot a_i \cdot z_i \cdot u_i$$

This product of two continuous decision variables makes the problem nonlinear.

---

## The McCormick Relaxation

The McCormick envelope replaces the exact bilinear constraint with four linear inequalities. 
Given bounds $z^L_i \leq z_i \leq z^U_i$ and $u^L_i \leq u_i \leq u^U_i$, where:

$$z^L_i = (V^{eff,min}_i)^{b_i}, \qquad z^U_i = (V^{eff,max}_i)^{b_i}$$

the four constraints that together form the envelope of $p_i = c \cdot a_i \cdot z_i \cdot u_i$ are:

$$\frac{p_i}{c \cdot a_i} \geq z^L_i \cdot u_i + u^L_i \cdot z_i - z^L_i \cdot u^L_i \tag{MC1}$$

$$\frac{p_i}{c \cdot a_i} \geq z^U_i \cdot u_i + u^U_i \cdot z_i - z^U_i \cdot u^U_i \tag{MC2}$$

$$\frac{p_i}{c \cdot a_i} \leq z^U_i \cdot u_i + u^L_i \cdot z_i - z^U_i \cdot u^L_i \tag{MC3}$$

$$\frac{p_i}{c \cdot a_i} \leq z^L_i \cdot u_i + u^U_i \cdot z_i - z^L_i \cdot u^U_i \tag{MC4}$$

MC1 and MC2 form the lower envelope (underestimates). MC3 and MC4 form the upper envelope
(overestimates). Together they form the **tightest possible linear relaxation** of the bilinear
product over the box $[z^L, z^U] \times [u^L, u^U]$.

The auxiliary variable $z_i$ is linked to $V_i$ via a Gurobi general power constraint:

$$z_i = V_i^{b_i} \tag{GP}$$

With `FuncNonlinear = 1`, Gurobi enforces this exactly using spatial branch-and-bound.

---

## Known Limitation: The Relaxation Gap

The McCormick envelope is a **relaxation**, not an exact formulation. For any point
$(z_i, u_i)$ strictly inside the box, the upper bound on $p_i$ from MC3 exceeds the true value:

$$\text{gap} = (z^U_i - z_i)(u_i - u^L_i) \geq 0$$

This means the solver can find solutions where $p^{relaxed}_i > p^{physical}_i$.
The physical power is always computed post-solve as:

$$p^{phys}_i = c \cdot a_i \cdot V_i^{b_i} \cdot u_i$$

and stored in `results.p`. The relaxation gap is quantified by the diagnostics in `mc_diagnostics.m`.

---

## Variable Layout

| Block | Variables | Size |
|---|---|---|
| 0 | $V_i$ — reservoir storage | $n$ |
| 1 | $p_i$ — power output | $n$ |
| 2 | $u_i$ — turbine release | $n$ |
| 3 | $s_i$ — spill | $n$ |
| 4 | $d_i$ — tracking error | $n$ |
| 5 | $z_i$ — head proxy $V_i^{b_i}$ | $n$ |

Total: $6n$ variables. All continuous. No binary variables.

---

## Objective Function

$$\min \quad -\sum_{i=1}^n p_i \;+\; 10^{-4} \sum_{i=1}^n s_i \;+\; \theta \sum_{i=1}^n \frac{d_i}{V^{max}_i - V^{min}_i}$$

---

## Constraint Summary

For each unit $i$ at time $t$:

| Label | Constraint | Type |
|---|---|---|
| C1 | $V_i + u_i + s_i = V^{prev}_i + \hat{q}_i$ | Mass balance |
| C2 | $u_i \geq u^{prev}_i + \Delta^{dn}_i$ | Ramp rate lower |
| C3 | $u_i \leq u^{prev}_i + \Delta^{up}_i$ | Ramp rate upper |
| MC1–MC4 | McCormick envelope of $p_i = c a_i z_i u_i$ | Power production |
| GP | $z_i = V_i^{b_i}$ | General power constraint |
| C5 | $V_i - d_i \leq V^{ref,hi}_i$ | Tracking upper |
| C6 | $-V_i - d_i \leq -V^{ref,lo}_i$ | Tracking lower |

For $t = 1$, C2 and C3 are omitted.

---

## Stochastic Bound Tightening

When `bounds = 'jcc-bon'`, the storage bounds are tightened using a Bonferroni correction.
The risk budget $\varepsilon$ is split equally across all $n$ units and both bound directions,
giving a safety margin of:

$$z = \Phi^{-1}\!\left(1 - \frac{\varepsilon}{2n}\right)$$

$$V^{eff,max}_i = V^{max}_i - z \cdot \hat{\sigma}_i, \qquad
V^{eff,min}_i = \max\!\left(0,\; V^{min}_i + z \cdot \hat{\sigma}_i\right)$$

The non-negativity guard on $V^{eff,min}$ is required because $z_i = V_i^{b_i}$ is undefined
for $V_i \leq 0$.

---

## Solver Configuration

```matlab
params.FuncNonlinear = 1;   % exact enforcement of z = V^b via spatial B&B
params.Threads       = 1;   % deterministic latency for real-time loop
params.TimeLimit     = 10;  % seconds — reduce to match control interval
params.Seed          = 1;   % reproducibility
```

---

## `sys` Struct Required Fields

| Field | Description |
|---|---|
| `.min_V`, `.max_V` | Physical storage bounds |
| `.a` | Head-storage coefficient ($a_i = h^{min}_i / (V^{min}_i)^{b_i}$ if not set) |
| `.b` | Head-storage exponent $\in (0, 1)$ |
| `.min_ut`, `.max_ut` | Turbine release bounds |
| `.RR_dn`, `.RR_up` | Ramp rate limits |
| `.F` | Installed power capacity |
| `.min_h`, `.max_h` | Head bounds (used to derive `.a` if absent) |

---

## Comparison with PWL-MILP

| Property | McCormick (`realtimeGurobiMC`) | PWL-MILP (`realtimeGurobiPWL`) |
|---|---|---|
| Problem class | MINLP (spatial B&B) | MILP (standard B&B) |
| Power overstatement | Yes — structural relaxation gap | No — approximation only |
| Binary variables | None | $nK$ |
| Accuracy control | Spatial B&B tolerance | Increase $K$ |
| `FuncNonlinear` | Required | Not needed |
| Use case | Benchmarking, gap analysis | Production real-time dispatch |

---

*For the production formulation see `realtimeGurobiPWL.m`.  
For the multi-period oracle see `oracleGurobi.m`.*

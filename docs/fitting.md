# Inflow forecast and conditional variance

The controller requires two quantities from the inflow record: a mean
forecast of local inflow and a conditional variance for its error. This
document gives the calculations that produce both, in the order they are
performed, and the relation between the estimated coefficients and the form
in which they are reported.

---

## Contents

1. [Inflow reconstruction](#1-inflow-reconstruction)
2. [Mean forecast](#2-mean-forecast)
3. [Conditional variance](#3-conditional-variance)
4. [Normalisation of the exogenous inputs](#4-normalisation-of-the-exogenous-inputs)
5. [Estimation](#5-estimation)
6. [Reported coefficients](#6-reported-coefficients)
7. [Use in the controller](#7-use-in-the-controller)

---

## 1. Inflow reconstruction

The USACE record gives hourly generation flow, spill and forebay elevation at
each project; it does not give inflow. Inflow is reconstructed from mass
balance.

Total release is the sum of both outlets,

```
O(i,t) = Flow-Gen(i,t) + Flow-Spill(i,t)
```

Spill is included throughout. Over a full water year it is non-zero in 40–72%
of hours and carries up to 54% of the water in May, so a generation-only
release both understates reconstructed inflow and removes most of the routed
upstream release from the exogenous inputs of section 4.

Storage is recovered from forebay elevation by inverting the project's affine
head relation and rescaling onto its usable range, giving `V(i,t)` in kaf.
Total inflow then follows from the hourly volume change,

```
q(i,t) = O(i,t) + ( V(i,t) - V(i,t-1) ) / k_af ,      k_af = 0.0826446
```

with `k_af` in kaf per kcfs-hour, so that the storage term is expressed as a
flow. Each reconstructed series is then shifted by a constant so that its
mean equals the mean of total release over the record, which closes mass over
the eight years and removes the accumulated bias of the elevation-to-storage
inversion.

The record is grouped into eight season years rolling over on 1 May, so that
each year contains one complete freshet. The recursion of section 3 resets at
every boundary.

---

## 2. Mean forecast

Local inflow at project `i` follows an AR(1), estimated by least squares:

```
mu(i,t) = a0(i) + a1(i) * q(i,t-1)
```

The routed upstream release enters the water balance of the dispatch program
explicitly, so the forecast carries no separate release term. McNary heads the
reach and draws its upstream inflow from the exogenous Priest Rapids record.

The forecast error is the residual of this regression,

```
e(i,t) = q(i,t) - mu(i,t)
```

and sections 3 to 5 model its variance. The coefficients `a0` and `a1` appear
in `results/tables/table_system_parameters.csv` as `ar0_kcfs` and `ar1`.

---

## 3. Conditional variance

The error variance is conditioned on its own recent magnitude and on the
upstream release decision:

```
sigma^2(i,t) = omega(i) + alpha(i) * e^2(i,t-1)
               + gu(i) * zu(j,t-tau) + gr(i) * zr(j,t-tau)
```

`j` is the upstream project and `tau` its travel lag. `zu` and `zr` are the
normalised upstream release level and the absolute hour-to-hour change in it,
defined in section 4. The two exogenous terms are what make the uncertainty
decision-dependent: the variance a project faces is a function of how the
project above it is being operated.

The specification is ARCH-X rather than GARCH-X. The persistence term
`beta(i) * sigma^2(i,t-1)` estimates at exactly zero for all four projects
and is therefore omitted; `beta` is carried as a zero column so that the
parameter set matches the general form.

At the first hour of each season year the recursion has no lagged residual,
and the variance is set to its unconditional value

```
sigmabar^2(i) = E[ e^2(i,t) ]
```

estimated as the sample mean of the squared residuals. The same quantity is
the constant that the decision-independent framework uses in place of the
whole expression above. It is stored as `sigma_DIU`.

---

## 4. Normalisation of the exogenous inputs

Both exogenous inputs are min–max normalised on their estimation range:

```
zu = ( u(j,t-tau) - u_min(j) ) / u_rng(j)
zr = ( |u(j,t-tau) - u(j,t-tau-1)| - r_min(j) ) / r_rng(j)
```

so each lies in `[0,1]` over the estimation sample. Outside that sample the
transform is applied and then clipped below at zero, so the model never
extrapolates beneath its estimation range. The four constants are stored as
`release_min`, `release_range`, `ramp_min` and `ramp_range`.

Two consequences matter for what follows. The coefficients `gu` and `gr` are
in kcfs² per unit of normalised input, so they are directly comparable with
each other and with `omega`. And the exogenous block attains its maximum
`gu + gr` at `zu = zr = 1`, that is at a full-scale upstream release
coinciding with a full-scale ramp.

McNary has no upstream project, so `zu` and `zr` are identically zero,
`gu = gr = 0` by construction, and the stored ranges are `(0,1)`.

---

## 5. Estimation

The variance parameters are estimated by Gaussian maximum likelihood on the
residuals, each project separately. With `theta = (omega, alpha, gu, gr)`,

```
-L(theta) = (1/2) * sum_t [ log(2*pi) + log sigma^2(i,t)
                            + e^2(i,t) / sigma^2(i,t) ]
```

summed over every hour except the first of each season year, where the
variance is the reset value of section 3 and carries no information about
`theta`. The parameters are bounded

```
omega > 0,      0 <= alpha <= 0.995,      gu >= 0,      gr >= 0
```

and `-L` is minimised by L-BFGS-B from two starting points, the better
solution being taken.

The bound on `alpha` keeps the recursion stationary. The non-negativity
bounds keep `sigma^2` positive for every admissible input, which the chance
constraint requires: a negative coefficient could drive the variance below
zero at an extreme of the normalised range and leave the safety margin
undefined.

The fit uses the whole eight-year record rather than a seasonal window:
70,080 hourly observations per project, with no train/test split. The full
year is necessary because the storage violations the constraint exists to
control concentrate in May and June, while a window drawn from the low-flow
season contains no spill, and would identify the exogenous coefficients on
the months in which the constraint never binds.

At John Day the release coefficient estimates at its lower bound, `gu = 0`,
and the decision dependence there is carried entirely by the ramp term.

The resulting coefficients are in `data/archx_params.csv`, one row per
project, and are read by the controller as given.

---

## 6. Reported coefficients

The coefficients are estimated in the units of the data: kcfs² for `omega`,
`gu` and `gr`, dimensionless for `alpha`. They are reported scaled by the
unconditional variance, with a single exogenous coefficient in place of two:

```
omega_reported = omega / sigmabar^2
alpha_reported = alpha
gamma_reported = ( gu + gr ) / sigmabar^2
```

The scaling makes every coefficient dimensionless and readable as a multiple
of the unconditional variance.

The aggregation is exact rather than a summary. Because `zu` and `zr` share
the `[0,1]` scale of section 4, `gu + gr` is the exogenous block evaluated at
`zu = zr = 1`, so `gamma_reported` is the largest multiple of `sigmabar^2`
that the upstream decision can contribute. The split between level and ramp
is retained in the stored coefficients, and the controller applies the two
separately.

| Project | `omega/sigmabar^2` | `alpha` | `(gu+gr)/sigmabar^2` |
|---|---|---|---|
| MCN | 0.632 | 0.441 | 0.00 |
| JDA | 0.447 | 0.239 | 5.17 |
| TDA | 0.100 | 0.061 | 10.73 |
| BON | 0.190 | 0.352 | 8.79 |

Two readings follow. `omega/sigmabar^2` is below one at every project, so the
baseline conditional variance sits beneath the unconditional one and the
remainder is supplied by the lagged residual and the upstream decision; that
is what distinguishes a conditional model from a constant. And the exogenous
range is wide: at The Dalles the decision terms alone can raise the variance
to roughly eleven times its unconditional value.

Both forms appear in `results/tables/table_archx_fits.csv`, the estimated
coefficients and the scaled ones as `omega_norm` and `gamma_norm`.

---

## 7. Use in the controller

`forecastError.m` evaluates section 3 one look-ahead step at a time, with the
exogenous block scaled by `Kappa`:

```
sigma^2 = omega(kappa) + alpha * e^2 + kappa * ( gu * zu + gr * zr )
```

`Kappa = 1` applies the estimated relationship and `Kappa = 0` gives the
decision-independent arm. The intercept is not `omega` itself but

```
omega(kappa) = omega + (1 - kappa) * ( gu * zu_bar + gr * zr_bar )
```

where `zu_bar` and `zr_bar` are the mean normalised inputs over the evaluated
year. As `kappa` falls to zero the exogenous block is replaced by its own
sample mean and folded into the intercept, so the two arms carry the same
average variance and differ only in whether that variance responds to the
decision. Without the adjustment the decision-independent arm would also be
the lower-variance arm, and a comparison between them would measure a
difference in level rather than in responsiveness.

Within the look-ahead the lagged squared residual is observed only at the
first step and is replaced thereafter by the preceding conditional variance,
the standard multi-step recursion. The returned standard deviation is clipped
to `[0.2, 3.0] * sigmabar`.

`docs/realtime.md` §4 covers how upstream ramps falling inside the look-ahead
are predicted without making the variance a function of the decision
variables.

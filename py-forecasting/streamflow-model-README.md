# Hourly Inflow and Decision-Dependent Error Model

## Scope

The notebook constructs hourly inflow series for the McNary (MCN), John Day (JDA), The Dalles (TDA), and Bonneville (BON) reservoirs, then estimates one-hour-ahead inflow models for JDA, TDA, and BON. The final uncertainty model separates predictable inflow dynamics, serial correlation in the remaining forecast errors, time-varying marginal variance, and contemporaneous cross-unit dependence.

The cascade is:

`MCN -> JDA -> TDA -> BON`

The analysis uses the September 7 through December 5 dry-season window. Training seasons are 2018-2022 and test seasons are 2023-2025.

## Inflow reconstruction

Hourly inflow is reconstructed from outflow and storage change using the reservoir mass balance:

`q_i,t^recon = u_i,t + (V_i,t - V_i,t-1) / K`

where `q_i,t^recon` is reconstructed inflow, `u_i,t` is total outflow, `V_i,t` is storage, and `K = 0.0826446 kaf/(kcfs-hour)` is the unit-conversion factor. Forebay elevation is first normalized and mapped to normalized storage using the inverse linear head-volume relationship.

A unit-specific mean bias is removed relative to the reported inflow series:

`q_i,t^bc = q_i,t^recon - mean(q_i,t^recon - q_i,t^reported)`

The active forecasting target in the notebook is a trailing six-hour rolling median of the bias-corrected reconstruction, denoted `q_i,t^avg`. This filter reduces high-frequency reconstruction noise while using only current and past observations.

## Travel-time selection

Candidate travel times are restricted using reservoir distance and plausible propagation speeds of 5-15 mph. For reach length `d`, the candidate lag set is:

`L = {ceil(d/15), ..., floor(d/5)}`

The selected travel time maximizes the correlation between lagged upstream outflow and downstream reconstructed inflow. The notebook reports 12 hours for MCN-JDA, 3 hours for JDA-TDA, and 7 hours for TDA-BON. The active forecasting specification instead uses the following empirically selected lag sets:

| Unit | Upstream unit | Inflow lags | Release lags |
|---|---|---|---|
| JDA | MCN | 1, 2, 3 | 9, 10, 11 |
| TDA | JDA | 1, 2, 3 | 1, 2, 3 |
| BON | TDA | 1 | 4 |

## Conditional mean model

For downstream unit `i` and upstream unit `j`, the one-hour-ahead conditional mean is:

`mu_i,t = b_i,0 + sum(phi_i,l q_i,t-l) + sum(theta_i,h u_j,t-h)`

Recent local inflows represent short-term persistence, while lagged upstream releases represent routed cascade effects. Each model is estimated by ridge regression. Predictors are standardized using training data only, and the ridge penalty is selected by cross-validation over `10^-4` to `10^4`. Ridge regularization is used because adjacent inflow and release lags are strongly collinear.

Testing is rolling one-step-ahead: observed inflow through hour `t-1` is used to predict hour `t`. This matches real-time dispatch, in which the latest realized inflow is available before the next decision. It also prevents artificial error accumulation from a 2,160-hour recursive forecast. Negative predictions are truncated to zero.

The resulting test coefficients of determination are 0.970 for JDA, 0.969 for TDA, and 0.983 for BON.

## Autoregressive error correction

Raw mean-model residuals are defined as:

`e_i,t = q_i,t - mu_i,t`

Residual diagnostics show remaining dependence at short lags and at the 24-hour operating cycle. A compact autoregressive correction is therefore estimated on training residuals:

`e_i,t = c_i + phi_i,1 e_i,t-1 + phi_i,2 e_i,t-2 + phi_i,24 e_i,t-24 + epsilon_i,t`

The innovation `epsilon_i,t` is the input to the variance model. The correction reduced test residual RMSE by approximately 9.0% for JDA, 6.3% for TDA, and 29.5% for BON. These lags are retained because they address observed residual structure while adding only three state variables, all known at the real-time decision point.

## Student-t GARCH-X variance model

The conditional innovation model is:

`epsilon_i,t = sigma_i,t z_i,t`

`z_i,t ~ standardized Student-t(nu_i)`

The decision-dependent conditional variance is:

`sigma_i,t^2 = omega_i + alpha_i epsilon_i,t-1^2 + beta_i sigma_i,t-1^2 + gamma_i,u u_tilde_j,t-h + gamma_i,r abs(Delta u_tilde_j,t-h)`

Here, `u_tilde` and `abs(Delta u_tilde)` are upstream release and absolute hourly ramp normalized to the training range. The lag `h` is the center of the active release-lag set: 10 hours for JDA, 2 hours for TDA, and 4 hours for BON.

Release level represents the operating regime, while absolute ramp represents the magnitude of an operational change. Absolute ramp is used because both upward and downward changes may increase routing uncertainty. Including both variables provided the best training fit and maintained strong out-of-sample likelihood performance.

Student-t innovations are used because the residual histograms and quantile plots show heavy tails, and Student-t GARCH produced substantially lower AIC than Gaussian GARCH. The estimated degrees of freedom are approximately 3.04 for JDA, 3.35 for TDA, and 2.64 for BON.

Parameters are estimated by maximum likelihood with `omega`, `alpha`, `beta`, `gamma_u`, and `gamma_r` constrained to be nonnegative. The stationarity safeguard is:

`alpha_i + beta_i <= 0.995`

The variance recursion is reset at the beginning of every seasonal year so that the end of one dry season is not treated as adjacent to the beginning of the next.

## Final fitted results

| Unit | phi_1 | phi_2 | phi_24 | alpha | beta | gamma_release | gamma_ramp | nu | Persistence |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| JDA | 0.0180 | -0.2009 | 0.3147 | 0.3303 | 0.6626 | 1.6326 | 22.6097 | 3.0435 | 0.9929 |
| TDA | 0.0008 | -0.1966 | 0.3495 | 0.3646 | 0.5078 | 7.1607 | 23.0376 | 3.3522 | 0.8724 |
| BON | 0.8454 | -0.3347 | 0.1912 | 0.5434 | 0.4516 | 0.8876 | 16.4993 | 2.6415 | 0.9950 |

Because release and ramp are normalized to the same 0-1 range, their variance coefficients are comparable within each unit. Ramp magnitude has the larger estimated effect for all three units. BON remains at the imposed persistence boundary and should be included in a persistence-cap sensitivity analysis.

Held-out test performance is:

| Unit | Mean NLL | 90% coverage | 95% coverage | 99% coverage |
|---|---:|---:|---:|---:|
| JDA | 2.6573 | 0.8868 | 0.9385 | 0.9925 |
| TDA | 2.7488 | 0.8897 | 0.9498 | 0.9939 |
| BON | 1.7686 | 0.8942 | 0.9344 | 0.9770 |

TDA is well calibrated. JDA is slightly narrow at the 90% and 95% levels. BON remains under-covered in the upper tail, despite the Student-t specification.

## Cross-unit covariance

The constant correlation matrix is estimated from standardized training innovations `z_i,t = epsilon_i,t / sigma_i,t`, not from raw residuals:

```text
          JDA       TDA       BON
JDA    1.0000    0.0400    0.0194
TDA    0.0400    1.0000    0.0785
BON    0.0194    0.0785    1.0000
```

At time `t`, the covariance matrix is reconstructed as:

`Sigma_t = D_t R D_t`

`D_t = diag(sigma_JDA,t, sigma_TDA,t, sigma_BON,t)`

The weak off-diagonal correlations indicate that most remaining dependence is represented through the unit-specific conditional means and decision-dependent variances rather than simultaneous shocks.

## Reproducibility notes

- The residual-model cells require both `train_results` and `test_results`, but the current notebook explicitly constructs only `test_results`. Training one-step predictions must be created with `train_years` before a clean restart can reproduce the GARCH-X fit.
- The active target is `*_inflow_avg`, the trailing six-hour median. To use the unsmoothed hourly reconstruction, replace it consistently with `*_inflow_recon_bc` in the travel-time and forecasting sections.
- The line that normalizes `inflow_data_norm` is commented out; the exported `inflow.csv` therefore contains unnormalized flow despite the variable name.
- The variable `V_kaf` is currently computed from normalized storage without multiplication by physical storage limits. Before interpreting reconstructed inflow in kcfs, normalized storage must be converted to physical kaf consistently with `K`.
- The travel-time table and active `reach_map` do not currently use identical lags. The final lag choice should be reconciled and documented before publication.
- Backward filling can use future information. Missing observations should preferably be filled within each seasonal year using a causal method.
- The name `fit_model` is used first for ridge regression and later for AR-GARCH-X estimation. Distinct names would make a clean notebook restart safer.

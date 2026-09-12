# Streamflow Forecasting Model

## Purpose

This workflow estimates hourly inflow at each modeled reservoir from two sources of information:

1. recent inflow at the same reservoir; and
2. current and lagged outflow from the next upstream reservoir.

The model is intended to represent water propagation through the cascade while retaining short-term persistence in local inflow. A separate ridge-regression model is trained for Bonneville (BON), The Dalles (TDA), and John Day (JDA). McNary is not predicted because the current dataset does not include an additional upstream release series for use as its routing input.

## Modeled reaches

| Predicted unit | Target inflow | Upstream release predictor |
|---|---|---|
| Bonneville (BON) | `bon_inflow_avg` | `tda_outflow` |
| The Dalles (TDA) | `tda_inflow_avg` | `jda_outflow` |
| John Day (JDA) | `jda_inflow_avg` | `mcn_outflow` |

The inflow and outflow tables are joined using an inner merge on `datetime`. Consequently, only timestamps present in both tables enter the modeling dataset. Rows are sorted chronologically, and missing values are currently filled using backward filling before model construction.

## Training and test periods

The data are divided by season rather than by randomly sampled observations:

- **Training seasons:** 2018–2022
- **Test seasons:** 2023–2025

A `season_year` label prevents lagged predictors from crossing seasonal boundaries. It is defined as

$$
\text{season\_year}_t = \text{year}(t)-\mathbb{1}\{\text{month}(t)\leq 4\}.
$$

Thus, observations from January through April are assigned to the season that began in the preceding calendar year. For the September–December study window, `season_year` is equal to the calendar year.

## Model structure

For each reservoir, the normalized inflow prediction is

$$
\widehat{q}^{\,*}_t
= \beta_0
+ \sum_{j=1}^{p}\beta_j q^{*}_{t-j}
+ \sum_{k=0}^{K}\gamma_k r^{*}_{t-k},
$$

where:

- $q_t$ is inflow at the predicted reservoir;
- $r_t$ is outflow from the next upstream reservoir;
- $p$ is the number of autoregressive inflow lags;
- $K$ is the maximum upstream-release lag;
- $q^*$ and $r^*$ denote normalized variables; and
- the $k=0$ term includes the upstream release at the current time step.

Because the data are hourly, each lag represents one hour. For example, $(p,K)=(6,12)$ uses the previous six hours of local inflow and upstream releases from the current hour through 12 hours earlier.

The tested lag structures are

```python
lag_pairs = [(3, 6), (6, 6), (6, 12)]
```

## Lagged-data construction

Lagged features are generated separately within each `season_year`:

- `q_lag_1, ..., q_lag_p` contain past inflows at the predicted reservoir;
- `r_lag_0, ..., r_lag_K` contain current and past releases from the upstream reservoir.

Rows without a complete lag history are removed. Therefore, the first $\max(p,K)$ observations of each season are unavailable for model fitting and evaluation.

## Normalization

All scaling parameters are estimated from the training subset only. Target inflows are normalized as

$$
q_t^*=\frac{q_t-q_{\min}}{q_{\max}-q_{\min}},
$$

and upstream releases are normalized as

$$
r_t^*=\frac{r_t-r_{\min}}{r_{\max}-r_{\min}}.
$$

The inflow and release minima and ranges are stored with the fitted model and reused without recalculation during testing. This prevents the test-period magnitude range from influencing the fitted coefficients or scaling parameters.

## Ridge-regression training

One ridge-regression model is fitted for every combination of reservoir and candidate lag structure. Ridge regression estimates the coefficient vector by minimizing

$$
\sum_t\left(q_t^*-\widehat{q}_t^*\right)^2
+\alpha\lVert\boldsymbol{\theta}\rVert_2^2,
$$

where $\boldsymbol{\theta}$ contains the fitted lag coefficients and $\alpha$ controls the amount of coefficient shrinkage. Penalization is useful here because adjacent hourly lags are strongly correlated.

`RidgeCV` selects $\alpha$ from 40 logarithmically spaced candidates:

```python
np.logspace(-4, 4, 40)
```

The notebook does not explicitly provide a cross-validation splitter to `RidgeCV`; therefore, alpha selection follows the default behavior of the installed scikit-learn version. The final model is then fitted using all training observations for that reservoir and lag structure.

## Recursive test procedure

Testing is performed independently for each test season. At the beginning of a season, the first $\max(p,K)$ observed inflows provide the warm-up history. Predictions then proceed chronologically:

1. Construct the predictor vector from the most recently available inflow estimates and observed upstream releases.
2. Normalize the predictors using the training-period scaling parameters.
3. Predict normalized inflow with the fitted ridge model.
4. Transform the prediction back to the original flow scale.
5. Replace the inflow value in the prediction history with the new prediction so it can be used at subsequent time steps.
6. Impose the physical lower bound $\widehat q_t\geq 0$.

The recursive update is

$$
\widehat q_t=\max\left\{0,\ q_{\min}+q_{\mathrm{range}}\widehat q_t^*\right\}.
$$

The first prediction uses observed inflow values for its entire autoregressive history. As the recursion advances, those observations are progressively replaced by model predictions. After $p$ predicted time steps, all local-inflow lag terms are recursive predictions. Upstream-release predictors remain observed values throughout the current test implementation.

## Performance metrics

For each reservoir and lag structure, the notebook calculates the prediction error as

$$
e_t=\widehat q_t-q_t.
$$

It reports:

- **RMSE:** $\sqrt{N^{-1}\sum_t e_t^2}$;
- **MAE:** $N^{-1}\sum_t |e_t|$;
- **bias:** $N^{-1}\sum_t e_t$, where positive bias indicates overprediction; and
- **correlation:** Pearson correlation between observed and predicted inflow.

Metrics are pooled across the selected test seasons. The lag structure with the lowest test-period RMSE is labeled as the best model for each reservoir.

## Forecast-error correlation

The final calculation measures dependence among the errors at BON, TDA, and JDA. For a selected lag pair, residuals are aligned by timestamp:

$$
\varepsilon_{i,t}=q_{i,t}-\widehat q_{i,t},
$$

and their Pearson correlation matrix is computed. This matrix summarizes contemporaneous dependence among forecast errors and can be used to parameterize the spatial uncertainty structure in the downstream optimization model.

The current notebook calculates this matrix using $(p,K)=(6,6)$ for all three reservoirs, rather than using the separately selected best lag structure for each reservoir.

## Interpretation and current limitations

The reported results should be interpreted with the following implementation details in mind:

1. **Lag selection uses the test period.** The 2023–2025 seasons are used both to select $(p,K)$ and to report performance. The resulting RMSE is therefore a model-selection score, not a fully independent estimate of out-of-sample performance. A stricter design would select $p$, $K$, and $\alpha$ using only 2018–2022, then evaluate the chosen specification once on 2023–2025.
2. **Backward filling can introduce future information.** `model_df.bfill()` replaces a missing value with a later observation and may also fill across seasonal boundaries. Missing-data treatment should be restricted by season and designed to avoid look-ahead leakage.
3. **Upstream releases are known during testing.** The test procedure uses observed $r_t$, including `r_lag_0`. This is appropriate only when the current upstream release is observable or supplied by the dispatch model at prediction time. It is not a fully autonomous multi-step forecast of both inflow and release.
4. **Warm-up observations are required.** Each seasonal simulation begins with observed target inflows. Performance excludes the warm-up interval.
5. **Only nonnegativity is enforced.** Predictions are clipped at zero but are not capped at a historical or physical maximum.
6. **One linear relationship is fitted per reach.** The model does not currently include seasonal interactions, nonlinear routing behavior, reservoir state, or time-varying travel time.

## Reproducibility summary

For each candidate $(p,K)$ and each modeled reservoir, the workflow is:

```text
merge inflow and outflow data
→ assign season_year
→ construct within-season lags
→ retain 2018–2022 training rows
→ estimate training-only normalization
→ select ridge penalty and fit coefficients
→ recursively predict 2023–2025
→ calculate RMSE, MAE, bias, and correlation
→ select the lowest-RMSE lag structure by reservoir
→ calculate cross-reservoir residual correlations
```

The ridge model is deterministic for fixed input data, preprocessing, lag choices, and software behavior; no random seed is required by the current implementation.

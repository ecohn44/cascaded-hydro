# Hourly Inflow Forecasting and Decision-Dependent Uncertainty

## 1. Purpose and model structure

This notebook constructs hourly inflow data for the Lower Columbia River cascade and estimates a real-time uncertainty model for John Day (JDA), The Dalles (TDA), and Bonneville (BON). McNary (MCN) supplies the upstream release predictor for JDA but is not itself forecast because no additional upstream reservoir is included.

The modeled cascade is

$$
\mathrm{MCN}\rightarrow\mathrm{JDA}\rightarrow\mathrm{TDA}\rightarrow\mathrm{BON}.
$$

The final model has four layers:

$$
\text{inflow reconstruction}\rightarrow\text{conditional mean}\rightarrow\text{AR correction}\rightarrow\text{GARCH-X variance}.
$$

This decomposition distinguishes predictable inflow dynamics from serially correlated forecast errors, time-varying marginal uncertainty, and cross-unit dependence.

## 2. Raw data and operating limits: cells 0-10

Cells 0-7 load hourly inflow, generation outflow, spill, and forebay elevation for MCN, JDA, TDA, and BON. Missing values are backward-filled and the unit data frames are collected in a common dictionary.

For observed outflow $u_{i,t}$, cells 8-10 estimate empirical release and ramp limits:

$$
\underline u_i=\min_t u_{i,t},\qquad \overline u_i=\max_t u_{i,t},
$$

$$
\Delta u_{i,t}=u_{i,t}-u_{i,t-1},
$$

$$
\underline r_i=\min_t\Delta u_{i,t},\qquad \overline r_i=\max_t\Delta u_{i,t}.
$$

These empirical values characterize the historical operating region and can inform release and ramp constraints in the dispatch model.

## 3. Hydrograph construction: cells 11-15

The four reported inflow series are aligned by timestamp and converted from kcfs to cubic meters per second using

$$
1\ \mathrm{kcfs}=28.316846592\ \mathrm{m^3/s}.
$$

The system-average inflow is

$$
\bar q_t=\frac{1}{4}\sum_{i\in\{\mathrm{MCN,JDA,TDA,BON}\}}q_{i,t}.
$$

Hourly observations are aggregated to daily means. For each calendar day $d$, the notebook computes the cross-year minimum, maximum, mean, and empirical quantiles $Q_{0.10}(d)$, $Q_{0.30}(d)$, $Q_{0.70}(d)$, and $Q_{0.90}(d)$. Leap days are removed so all years share a common 365-day plotting axis. The resulting hydrograph shows both interannual trajectories and the historical operating envelope.

## 4. Dry-season definition and ranking: cells 16-23

Cells 16-21 use The Dalles spill to identify extended low-spill periods. For the 14-day window $W=336$ hours, the implemented rolling statistic is

$$
S_t=\sum_{k=t-W+1}^{t}s_k,
$$

where $s_k$ is hourly spill. The code classifies a period as dry when $S_t<168$. Contiguous dry blocks shorter than 40 days are discarded.

The fixed modeling season used later is September 7 through December 5. For each year $y$, cells 22-23 rank this 90-day period using

$$
\bar q_y=\frac{1}{T_y}\sum_{t\in\mathcal T_y}\bar q_t,
$$

where $\mathcal T_y$ is the hourly dry-season index. Years are ordered from smallest to largest $\bar q_y$ to identify the driest historical operating conditions.

## 5. Hourly inflow and storage data sets: cells 24-47

For years 2018-2025, the expected dry-season index is

$$
\mathcal T=\bigcup_{y=2018}^{2025}\{\text{Sep. 7 00:00 through Dec. 5 23:00 of year }y\}.
$$

Duplicate timestamps are averaged, missing timestamps are inserted, and missing values are backward-filled. With 90 days per year, the expected sample size is

$$
8\times90\times24=17{,}280\ \text{hours}.
$$

The notebook defines a common inflow min-max transformation

$$
q_{i,t}^{\mathrm{norm}}=\frac{q_{i,t}-q_{\min}}{q_{\max}-q_{\min}},
$$

where $q_{\min}$ and $q_{\max}$ are calculated across all units and dry-season observations. This preserves relative flow magnitudes across the cascade. In the current notebook, however, the line applying this transformation is commented out; the exported `inflow.csv` therefore contains unnormalized flow.

Forebay elevation is normalized separately for each reservoir:

$$
H_{i,t}^{\mathrm{norm}}=\frac{H_{i,t}-\min_tH_{i,t}}{\max_tH_{i,t}-\min_tH_{i,t}}.
$$

Unit-specific normalization is used because each reservoir has a different forebay operating range. The resulting storage-reference data are exported to `norm-soc.csv`.

## 6. Inflow reconstruction: cells 48-64

The notebook reconstructs hourly inflow from outflow and storage change. The fitted normalized head-volume relationship is

$$
H_{i,t}^{\mathrm{norm}}=a_iV_{i,t}^{\mathrm{norm}}+b_i,
$$

which is inverted as

$$
V_{i,t}^{\mathrm{norm}}=\mathrm{clip}(\frac{H_{i,t}^{\mathrm{norm}}-b_i}{a_i},0,1).
$$

For hourly outflow $u_{i,t}$ and physical storage $V_{i,t}$ in kaf, mass balance gives

$$
q_{i,t}^{\mathrm{recon}}=u_{i,t}+\frac{V_{i,t}-V_{i,t-1}}{K},
$$

where

$$
K=0.0826446\ \frac{\mathrm{kaf}}{\mathrm{kcfs\cdot hour}}.
$$

A unit-specific additive bias is estimated relative to reported inflow:

$$
b_i^{q}=\frac{1}{N_i}\sum_t(q_{i,t}^{\mathrm{recon}}-q_{i,t}^{\mathrm{reported}}),
$$

$$
q_{i,t}^{\mathrm{bc}}=q_{i,t}^{\mathrm{recon}}-b_i^{q}.
$$

The active forecasting target is a trailing six-hour median:

$$
q_{i,t}^{\mathrm{avg}}=\mathrm{median}\{q_{i,t-k}^{\mathrm{bc}}:k=0,\ldots,5\}.
$$

The trailing window is causal because it contains no future observations. Relative to reported inflow, the six-hour series has RMSE values of 13.89, 19.49, 18.80, and 20.47 kcfs for BON, TDA, JDA, and MCN, respectively.

## 7. Travel-time estimation: cells 65-66

For reach distance $d_{ji}$ and plausible velocity interval $[v_{\min},v_{\max}]=[5,15]$ mph, the candidate travel-time set is

$$
\mathcal L_{ji}=\{\lceil d_{ji}/v_{\max}\rceil,\ldots,\lfloor d_{ji}/v_{\min}\rfloor\}.
$$

The selected lag maximizes the sample correlation between upstream outflow and downstream reconstructed inflow:

$$
h_{ji}^{*}=\underset{h\in\mathcal L_{ji}}{\mathrm{arg\,max}}\ \mathrm{Corr}(u_{j,t-h},q_{i,t}^{\mathrm{avg}}).
$$

The notebook reports 12 hours for MCN-JDA, 3 hours for JDA-TDA, and 7 hours for TDA-BON, corresponding to correlations of 0.679, 0.812, and 0.735.

## 8. One-hour-ahead conditional mean: cells 67-85

Training seasons are 2018-2022 and test seasons are 2023-2025. Lagged predictors are generated separately within each seasonal year, preventing the end of one dry season from being linked to the next.

For downstream reservoir $i$ and upstream reservoir $j$, the conditional mean is

$$
\mu_{i,t}=\theta_{i,0}+\sum_{\ell\in\mathcal P_i}\phi_{i,\ell}q_{i,t-\ell}^{\mathrm{avg}}+\sum_{h\in\mathcal L_i}\psi_{i,h}u_{j,t-h}.
$$

The active lag sets are:

| Unit | Upstream unit | $\mathcal P_i$ | $\mathcal L_i$ |
|---|---|---|---|
| JDA | MCN | 1, 2, 3 | 9, 10, 11 |
| TDA | JDA | 1, 2, 3 | 1, 2, 3 |
| BON | TDA | 1 | 4 |

The model is estimated by ridge regression:

$$
\widehat{\boldsymbol\theta}_i=\underset{\boldsymbol\theta_i}{\mathrm{arg\,min}}\ \{\lVert\mathbf y_i-\mathbf X_i\boldsymbol\theta_i\rVert_2^2+\lambda_i\lVert\boldsymbol\theta_i\rVert_2^2\}.
$$

Predictors are standardized using training data only, and `RidgeCV` selects $\lambda_i$ from $10^{-4}$ through $10^4$. Ridge regularization stabilizes estimates when adjacent inflow and release lags are strongly collinear.

Testing is rolling one-step-ahead:

$$
\widehat q_{i,t}=\max\{0,\widehat\mu_{i,t}(q_{i,t-1}^{\mathrm{obs}},q_{i,t-2}^{\mathrm{obs}},\ldots)\}.
$$

Observed lagged inflow is appropriate because the latest realization is available before each real-time dispatch decision. It also avoids artificial error accumulation from recursively predicting an entire 2,160-hour season. The test coefficients of determination are 0.970 for JDA, 0.969 for TDA, and 0.983 for BON, where

$$
R_i^2=1-\frac{\sum_t(q_{i,t}-\widehat q_{i,t})^2}{\sum_t(q_{i,t}-\bar q_i)^2}.
$$

Cells 70 and 78 evaluate residual variance, marginal distribution, autocorrelation, RMSE, MAE, bias, and correlation. The raw residual is

$$
e_{i,t}=q_{i,t}-\widehat q_{i,t}.
$$

## 9. Autoregressive residual correction: cells 86-88

The raw residual ACF retains dependence at lags 1, 2, and 24. The notebook therefore estimates

$$
e_{i,t}=c_i+\varphi_{i,1}e_{i,t-1}+\varphi_{i,2}e_{i,t-2}+\varphi_{i,24}e_{i,t-24}+\varepsilon_{i,t}.
$$

The corrected one-hour forecast is

$$
\widehat q_{i,t}^{\mathrm{corr}}=\widehat q_{i,t}+\widehat e_{i,t}^{\mathrm{AR}},
$$

and the remaining innovation is

$$
\varepsilon_{i,t}=q_{i,t}-\widehat q_{i,t}^{\mathrm{corr}}.
$$

The lag-24 term represents the observed daily operating cycle, while lags 1 and 2 represent short-term residual dynamics. All three lagged errors are known when the next hourly forecast is issued. The correction reduced test RMSE by 9.0% for JDA, 6.3% for TDA, and 29.5% for BON.

## 10. Student-t GARCH-X model: cells 87-90

The corrected innovation is written as

$$
\varepsilon_{i,t}=\sigma_{i,t}z_{i,t},
$$

where $z_{i,t}$ follows a standardized Student-t distribution with $\nu_i>2$ and unit variance. The decision-dependent conditional variance is

$$
\sigma_{i,t}^2=\omega_i+\alpha_i\varepsilon_{i,t-1}^2+\beta_i\sigma_{i,t-1}^2+\gamma_{i,u}\widetilde u_{j,t-h_i}+\gamma_{i,r}|\Delta\widetilde u_{j,t-h_i}|.
$$

The exogenous inputs are normalized using training ranges:

$$
\widetilde u_{j,t}=\frac{u_{j,t}-u_j^{\min}}{u_j^{\max}-u_j^{\min}},
$$

$$
|\Delta\widetilde u_{j,t}|=\frac{|u_{j,t}-u_{j,t-1}|-r_j^{\min}}{r_j^{\max}-r_j^{\min}}.
$$

The GARCH-X travel lag is the center of the active release-lag set: $h_{\mathrm{JDA}}=10$, $h_{\mathrm{TDA}}=2$, and $h_{\mathrm{BON}}=4$ hours. Release level identifies the upstream operating regime. Absolute ramp measures the magnitude of an operational change and treats sharp upward and downward changes symmetrically. Both predictors are retained because they improve the training likelihood and preserve strong test performance.

Student-t innovations are selected because the residual distributions are sharply peaked and heavy-tailed, and Student-t GARCH produced substantially smaller AIC than Gaussian GARCH. For standardized Student-t innovations, the implemented log-likelihood contribution is

$$
\ell_{i,t}=\log\Gamma(\frac{\nu_i+1}{2})-\log\Gamma(\frac{\nu_i}{2})-\frac{1}{2}\log[\pi(\nu_i-2)]-\frac{1}{2}\log\sigma_{i,t}^2-\frac{\nu_i+1}{2}\log(1+\frac{\varepsilon_{i,t}^2}{(\nu_i-2)\sigma_{i,t}^2}).
$$

The parameters maximize $\sum_t\ell_{i,t}$ subject to

$$
\omega_i>0,\qquad \alpha_i,\beta_i,\gamma_{i,u},\gamma_{i,r}\geq0,
$$

$$
\alpha_i+\beta_i\leq0.995,\qquad \nu_i\geq2.05.
$$

The persistence measure is $\alpha_i+\beta_i$. Values near one indicate slowly decaying volatility shocks. The variance recursion is reset at every seasonal boundary, with initial variance equal to the training innovation variance, so missing months do not enter the hourly recursion.

## 11. Out-of-sample evaluation: cell 87

Test filtering uses the fitted training parameters and the previous realized innovation. The latter is available in one-step real-time operation. No test observation is used to estimate the AR coefficients, normalization ranges, or GARCH-X parameters.

For nominal two-sided coverage $p$, the standardized Student-t critical value is

$$
c_{p,\nu}=\sqrt{\frac{\nu-2}{\nu}}\,T_{\nu}^{-1}(\frac{1+p}{2}).
$$

The interval is

$$
\widehat q_{i,t}^{\mathrm{corr}}\pm c_{p,\nu_i}\sigma_{i,t},
$$

and empirical coverage is

$$
\widehat C_{i,p}=\frac{1}{N_i}\sum_t\mathbf 1\{|\varepsilon_{i,t}|\leq c_{p,\nu_i}\sigma_{i,t}\}.
$$

The reported mean negative log-likelihood is

$$
\mathrm{NLL}_i=-\frac{1}{N_i}\sum_t\ell_{i,t}.
$$

| Unit | Test NLL | 90% coverage | 95% coverage | 99% coverage |
|---|---:|---:|---:|---:|
| JDA | 2.6573 | 0.8868 | 0.9385 | 0.9925 |
| TDA | 2.7488 | 0.8897 | 0.9498 | 0.9939 |
| BON | 1.7686 | 0.8942 | 0.9344 | 0.9770 |

TDA is well calibrated. JDA is slightly under-covered at 90% and 95%. BON remains under-covered in the upper tail.

## 12. Estimated parameters

| Unit | $\varphi_1$ | $\varphi_2$ | $\varphi_{24}$ | $\alpha$ | $\beta$ | $\gamma_u$ | $\gamma_r$ | $\nu$ | Persistence |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| JDA | 0.0180 | -0.2009 | 0.3147 | 0.3303 | 0.6626 | 1.6326 | 22.6097 | 3.0435 | 0.9929 |
| TDA | 0.0008 | -0.1966 | 0.3495 | 0.3646 | 0.5078 | 7.1607 | 23.0376 | 3.3522 | 0.8724 |
| BON | 0.8454 | -0.3347 | 0.1912 | 0.5434 | 0.4516 | 0.8876 | 16.4993 | 2.6415 | 0.9950 |

Because release and ramp are both normalized to $[0,1]$, their coefficients are comparable within each unit. Ramp magnitude has the larger estimated variance effect for all three units. BON reaches the imposed persistence limit and should be subjected to persistence-cap sensitivity analysis.

## 13. Diagnostics and DDU figures: cells 88-90

The standardized innovation is

$$
z_{i,t}=\frac{\varepsilon_{i,t}}{\sigma_{i,t}}.
$$

The ACF of $z_{i,t}$ diagnoses remaining conditional-mean dependence, while the ACF of $z_{i,t}^2$ diagnoses remaining conditional-variance dependence. The squared-innovation ACF is close to zero after GARCH-X filtering, indicating that the model captures most variance clustering. Some low-amplitude periodic dependence remains in $z_{i,t}$.

Cell 89 compares the DDU interval $\pm c_{0.95,\nu_i}\sigma_{i,t}$ with the constant DIU interval $\pm c_{0.95,\nu_i}s_i$, where $s_i$ is the training innovation standard deviation. Cell 90 evaluates the fitted response surface over release and ramp while holding the previous squared innovation and variance at their median values. These figures show both the temporal adaptation and the decision-dependent mechanism.

## 14. Cross-unit correlation and covariance: cells 91-92

The constant correlation matrix is estimated from aligned standardized training innovations:

$$
R=\mathrm{Corr}(\mathbf z_t)=
\begin{bmatrix}
1 & 0.0400 & 0.0194\\
0.0400 & 1 & 0.0785\\
0.0194 & 0.0785 & 1
\end{bmatrix},
$$

with unit order $[\mathrm{JDA},\mathrm{TDA},\mathrm{BON}]$. Define

$$
D_t=\mathrm{diag}(\sigma_{\mathrm{JDA},t},\sigma_{\mathrm{TDA},t},\sigma_{\mathrm{BON},t}).
$$

The time-varying covariance matrix is

$$
\Sigma_t=D_tRD_t.
$$

The small off-diagonal correlations indicate that most modeled dependence enters through unit-specific conditional means and decision-dependent variances rather than strong contemporaneous shocks.

## 15. Reproducibility issues to resolve

- The residual-model cells require `train_results`, but the current notebook explicitly constructs only `test_results`. One-step predictions for `train_years` must be created before a clean restart can reproduce the GARCH-X estimates.
- The active target is the six-hour median `*_inflow_avg`. To use unsmoothed hourly inflow, replace it consistently with `*_inflow_recon_bc` in travel-time estimation, `reach_map`, and `model_df`.
- The code labels normalized storage as `V_kaf` without multiplying by physical storage limits. A physical conversion is required before combining storage change with `K` in kaf/(kcfs-hour).
- The reported correlation-optimal travel times $(12,3,7)$ hours do not match all active release lags in `reach_map`. The final specification should reconcile this difference.
- The inflow min-max normalization line is commented out, so `inflow.csv` is not normalized despite the variable name.
- Backward filling can use future observations. Missing data should preferably be filled causally and separately within each seasonal year.
- `fit_model` is used for both ridge estimation and AR-GARCH-X estimation. Distinct function names are recommended for a clean top-to-bottom run.

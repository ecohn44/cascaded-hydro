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

Relative to reported inflow, the six-hour series has RMSE values of 13.89, 19.49, 18.80, and 20.47 kcfs for BON, TDA, JDA, and MCN, respectively.

## 7. Travel-time estimation: cells 65-66

For reach distance $d_{ji}$ and plausible velocity interval $[v_{\min},v_{\max}]=[5,15]$ mph, the candidate travel-time set is

$$
\mathcal L_{ji}=\{\lceil d_{ji}/v_{\max}\rceil,\ldots,\lfloor d_{ji}/v_{\min}\rfloor\}.
$$

The selected lag maximizes the sample correlation between upstream outflow and downstream reconstructed inflow:

$$
h_{ji}^{*}=\underset{h\in\mathcal L_{ji}}{\mathrm{arg\,max}}\ \mathrm{Corr}(u_{j,t-h},q_{i,t}^{\mathrm{avg}}).
$$

The notebook reports 10 hours for MCN-JDA, 2 hours for JDA-TDA, and 4 hours for TDA-BON, corresponding to correlations of 0.679, 0.812, and 0.735.

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
\hat{\theta}_i=\underset{\theta_i}{\mathrm{arg\,min}}\ \{|\mathbf y_i-\mathbf X_i\theta_i|_2^2+\lambda_i|\theta_i|_2^2\}.
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


## 9. GARCH-X model: cells 87-90

The model assumes 

$$
e_{i,t}\mid\mathcal F_{t-1}\sim\mathcal N(0,\sigma_{i,t}^2).
$$

and models the  decision-dependent conditional variance as

$$
\sigma_{i,t}^2=\omega_i+\alpha_i e_{i,t-1}^2+\beta_i\sigma_{i,t-1}^2+\gamma_{i,u}\widetilde u_{j,t-h_i}+\gamma_{i,r}|\Delta\widetilde u_{j,t-h_i}|.
$$

The exogenous inputs are normalized using training ranges:

$$
\widetilde u_{j,t}=\frac{u_{j,t}-u_j^{\min}}{u_j^{\max}-u_j^{\min}},
$$

$$
|\Delta\widetilde u_{j,t}|=\frac{|u_{j,t}-u_{j,t-1}|-r_j^{\min}}{r_j^{\max}-r_j^{\min}}.
$$

The GARCH-X travel lag is the center of the active release-lag set: $h_{\mathrm{JDA}}=10$, $h_{\mathrm{TDA}}=2$, and $h_{\mathrm{BON}}=4$ hours. Release level identifies the upstream operating regime. Absolute ramp measures the magnitude of an operational change and treats sharp upward and downward changes symmetrically. Both predictors are retained because they improve the training likelihood and preserve strong test performance.

For a candidate parameter vector

$$ 
\theta_i=
(\omega_i,\alpha_i,\beta_i,\gamma_{i,u},\gamma_{i,r}),
$$

the function calculates $\sigma_{i,t}^2$ recursively at every training hour. The resulting Gaussian log-likelihood contribution is

$$
\ell_{i,t}(\theta_i) = -\frac{1}{2}
\left[
\log(2\pi)
+\log\left(\sigma_{i,t}^2(\theta_i)\right)
+\frac{e_{i,t}^2}{\sigma_{i,t}^2(\theta_i)}
\right].
$$

Maximum likelihood estimation selects the parameters that maximize the sum of these contributions:

$$
\arg \max_{\theta_i} \sum_t\ell_{i,t}(\theta_i).
$$

The parameters maximize $\sum_t\ell_{i,t}$ subject to

$$
\omega_i>0,\qquad \alpha_i,\beta_i,\gamma_{i,u},\gamma_{i,r}\geq0,
$$

$$
\alpha_i+\beta_i\leq0.995.
$$

The persistence measure is $\alpha_i+\beta_i$. Values near one indicate slowly decaying volatility shocks. The variance recursion is reset at every seasonal boundary, with initial variance equal to the training innovation variance, so missing months do not enter the hourly recursion.

## 10. Out-of-sample evaluation: cell 87

Test filtering uses the GARCH-X parameters and normalization ranges estimated from the training data. At each test hour, the variance recursion uses the previous realized forecast residual, which is available in one-step-ahead real-time operation:

$$
e_{i,t-1}=q_{i,t-1}-\widehat q_{i,t-1}.
$$

At the beginning of each seasonal year, the recursion is reset using the training residual variance. Test observations are used only to update the one-step-ahead variance and evaluate performance; they are not used to re-estimate the normalization ranges or GARCH-X parameters.

Under the conditional Gaussian assumption, the critical value for a nominal two-sided coverage level \(p\) is

$$
c_p=\Phi^{-1}\left(\frac{1+p}{2}\right),
$$

where $\Phi^{-1}$ is the standard-normal quantile function. The prediction interval is

$$
\widehat q_{i,t}\pm c_p\sigma_{i,t}.
$$

Empirical coverage is calculated as

$$
\widehat C_{i,p}
=
\frac{1}{N_i}
\sum_t
\mathbf 1
\left\{
|e_{i,t}|
\leq
c_p\sigma_{i,t}
\right\}.
$$


| Unit | Test NLL | 90% coverage | 95% coverage | 99% coverage |
|---|---:|---:|---:|---:|
|JDA	|3.3482 | 0.9017|0.9323 |0.9691|
|TDA    |3.3685	|0.9051	|0.9389	|0.9753|
|BON	|2.8223	|0.9230	|0.9410	|0.9632|

The 90% Gaussian intervals are close to or slightly above their nominal coverage. All three units are under-covered at the 95% and 99% levels, indicating that the Gaussian model produces intervals that are too narrow in the tails. BON has the largest 90% coverage but the greatest under-coverage at the 99% level.

## 11. Estimated parameters

| Unit | $\omega$ | $\alpha$ | $\gamma_u$ | $\gamma_r$ | $\sigma_{\mathrm{DIU}}$ |
|---|---:|---:|---:|---:|---:|
| JDA | 26.5990 | 0.1605 | 47.3650 | 30.6786 | 7.03 |
| TDA | 25.2339 | 0.1576 | 62.8574 | 27.9373 | 7.28 |
| BON | 12.0440 | 0.2788 | 1.0294 | 30.6790 | 4.26 |

Because release and ramp are both normalized to $[0,1]$, their coefficients are comparable within each unit. Ramp magnitude has the larger estimated variance effect for all three units. 


## 12. Cross-unit correlation

The constant correlation matrix is estimated from aligned standardized training innovations:

$$
R =
\left[
\begin{array}{ccc}
1 & 0.09 & 0.068 \\
0.09 & 1 & 0.014 \\
0.068 & 0.046 & 1
\end{array}
\right]
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

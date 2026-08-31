# Kernel SOC Reference Learning

This folder contains the minimum implementation of the Kernel-1 reference-learning method in Qi et al. (2026), adapted from net load and hydrogen storage to normalized inflow and multi-unit hydropower SOC.

## File

- `kernel_soc_reference.m`: computes the next SOC reference and the weight assigned to every historical scenario.

There is no separate fitted model. Kernel regression is a lazy learner: the historical inflow and oracle SOC trajectories are the training data. Offline training only selects the window length `window` and Gaussian bandwidth `sigma`, which can be tuned later.

## Data

Stack each historical scenario along the third array dimension:

```matlab
size(inflow_norm)   % T x n x S
size(soc_ref)       % T x n x S
```

where:

- `T` is the number of time rows.
- `n` is the number of hydropower units.
- `S` is the number of historical scenarios.
- `soc_ref(1,:,s)` is the initial SOC for scenario `s`.
- `inflow_norm(t,:,s)` is the normalized inflow observed before the transition from SOC row `t` to row `t+1`.

At real-time step `t`, the online inputs are:

```matlab
size(inflow_history)   % t x n
size(soc_history)      % t x n
```

`soc_history` should contain measured reservoir SOC, including the current value. It should not contain future oracle references.

## Methods

Let the rolling-window index set be

$$
\mathcal W_t=\{\max(1,t-W+1),\ldots,t\},
$$

where $W$ is `window`. Vectorize and concatenate the observed inflow and SOC histories:

$$
\mathbf x_t = [ \mathbf \tilde{q}^{RT}, \mathbf v^{RT}  ]_{\mathcal{W_t}}
$$
 
For every historical scenario $s\in\{1,\ldots,S\}$, construct the corresponding training vector

$$
\mathbf y_{s,t}= \left[ \mathbf \tilde{q}^{(s)}, \mathbf v^{*,(s)}  \right]_{\mathcal W_t}
$$

Here, $\tilde{\mathbf q}$ is normalized inflow, $\mathbf v$ is measured online SOC, and $\mathbf v^{*,(s)}$ is the hindsight-optimal SOC trajectory for scenario $s$.

The feature-vector dimension is

$$
d_t=2n|\mathcal W_t|.
$$

The Gaussian similarity between the online history and scenario $s$ is

$$
K_{s,t}=\exp\left(
-\frac{\|\mathbf x_t-\mathbf y_{s,t}\|_2^2}
{2\sigma^2d_t}
\right),
$$

where $\sigma>0$ is the bandwidth. Normalize these similarities to obtain scenario weights:

$$
\omega_{s,t}=\frac{K_{s,t}}{\sum_{r=1}^{S}K_{r,t}},
\qquad
\sum_{s=1}^{S}\omega_{s,t}=1.
$$

The next multi-unit SOC reference is the weighted average of the historical oracle references:

$$
\hat{\mathbf v}_{t+1}
=\sum_{s=1}^{S}\omega_{s,t}\mathbf v_{t+1}^{*,(s)}.
$$

This is Qi et al.'s Kernel-1 construction with normalized inflow replacing normalized net load. The paper uses observations through $t-1$ to predict $\hat{\mathbf v}_t$; the MATLAB function uses the equivalent shifted convention of observations through row `t` to predict `soc_ref(t+1,:)`.

## Real-Time Use

```matlab
window = 168;   % Placeholder
sigma  = 0.10;  % Placeholder

[soc_target, weights] = kernel_soc_reference( ...
    inflow_norm, ...
    soc_ref, ...
    inflow_realized(1:t, :), ...
    soc_measured(1:t, :), ...
    window, ...
    sigma);
```

- `soc_target` is a `1 x n` SOC reference for time `t+1`.
- `weights` is an `S x 1` vector showing the contribution of each historical scenario.

Pass `soc_target` to the real-time dispatch reference-tracking term. After dispatch, append the newly measured inflow and SOC before computing the following reference.

## Tune `window` and `sigma`

Use leave-one-scenario-out cross-validation so the scenario being predicted is never included in its own training data:

```matlab
W_grid = [24 72 168 336];
sigma_grid = 0.01:0.01:0.10;

[best_W, best_sigma, mean_rmse, fold_rmse] = ...
    tune_kernel_parameters(inflow_norm, soc_ref, W_grid, sigma_grid);
```

For each candidate pair and each held-out scenario, the tuning function:

1. Trains on the other `S-1` scenarios.
2. Initializes the predicted trajectory with the held-out scenario's initial SOC.
3. Generates rows `2:T` recursively using observed inflow and previously predicted SOC.
4. Computes RMSE against the held-out oracle SOC.
5. Selects the pair with the lowest RMSE averaged across all held-out scenarios.

`mean_rmse(i,j)` is the average validation RMSE for `W_grid(i)` and `sigma_grid(j)`. `fold_rmse(i,j,s)` is the RMSE for held-out scenario `s` and should be inspected to ensure one year is not dominating the average.

## Historical Test Before Real-Time Dispatch

Reserve one complete scenario as the final test. For example, if scenario 8 is 2025:

```matlab
final_test = 8;
train = setdiff(1:size(inflow_norm, 3), final_test);

[best_W, best_sigma] = tune_kernel_parameters( ...
    inflow_norm(:, :, train), soc_ref(:, :, train), ...
    W_grid, sigma_grid);

[T, n, ~] = size(inflow_norm);
soc_test = zeros(T, n);
soc_test(1, :) = soc_ref(1, :, final_test);

for t = 1:T-1
    soc_test(t + 1, :) = kernel_soc_reference( ...
        inflow_norm(:, :, train), soc_ref(:, :, train), ...
        inflow_norm(1:t, :, final_test), soc_test(1:t, :), ...
        best_W, best_sigma);
end

test_error = soc_test(2:end, :) - soc_ref(2:end, :, final_test);
test_rmse = sqrt(mean(test_error(:).^2));
```

First compare `soc_test` with the held-out oracle trajectory visually and calculate both total RMSE and RMSE for each unit. The kernel result should also be compared with the simple historical-mean reference; otherwise the added kernel weighting has not demonstrated value.

Next run the real-time dispatch model on each held-out historical year. Pass the SOC produced by the dispatch simulation, rather than the held-out oracle SOC, into `kernel_soc_reference`. Record reference RMSE, terminal SOC error, SOC and flow constraint violations, generation, spill, and objective value. Do not feed the held-out oracle SOC history into the online function because that would make the backtest optimistic.

## Numerical Fallback

If all Gaussian similarities underflow to zero, the function assigns equal weights $1/S$. This prevents an undefined reference and reduces the prediction to the historical mean for that step.

## Reference

N. Qi, Y. Baker, and B. Xu, "Online dispatch of multi-duration energy storage in low-carbon microgrid with explainable reference learning," Applied Energy, vol. 412, 2026, Art. no. 127725.

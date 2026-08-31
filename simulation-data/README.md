# Cascaded Hydropower Data

This directory contains the raw and processed data used for the Bonneville, The Dalles, John Day, and McNary hydropower cascade. The raw flow and forebay-elevation time series were obtained from the U.S. Bureau of Reclamation (USBR). The processing workflow is implemented in `hydrograph.ipynb`.

## Reservoirs and file names

| Reservoir | USBR code | Raw flow file | Raw forebay file |
|---|---:|---|---|
| Bonneville | `BON` | `bonneville.csv` | `bonneville-soc.csv` |
| The Dalles | `TDA` | `dalles.csv` | `dalles-soc.csv` |
| John Day | `JDA` | `johnday.csv` | `johnday-soc.csv` |
| McNary | `MCN` | `mcnary.csv` | `mcnary-soc.csv` |

## Raw data

Each raw flow file contains a `Date Time` column and the following USBR series. Replace `XXX` with the reservoir code listed above.

| USBR column | Description | Units | Use in notebook |
|---|---|---:|---|
| `XXX.Flow-Gen.Ave.1Hour.1Hour.CBT-REV [kcfs]` | Hourly average generation flow | kcfs | Renamed `outflow (kcfs)` |
| `XXX.Flow-In.Inst.~6Hours.0.RFC-FCST [kcfs]` | Instantaneous inflow forecast, updated at approximately 6-hour intervals | kcfs | Renamed `inflow (kcfs)` |
| `XXX.Flow-Spill.Ave.1Hour.1Hour.CBT-REV [kcfs]` | Hourly average spill flow | kcfs | Renamed `spill (kcfs)` |

Each `*-soc.csv` file contains historical forebay elevation:

```text
XXX.Elev-Forebay.Inst.1Hour.0.CBT-REV [ft]
```

The notebook renames this field `soc_ft`. Forebay elevation is the historical measurement used to construct the model's scaled state-of-charge (SOC) trajectory; it is not a direct measurement of energy storage.

## Processed files

| File | Columns | Coverage | Description |
|---|---|---|---|
| `inflow.csv` | `datetime`, `bon_inflow`, `tda_inflow`, `jda_inflow`, `mcn_inflow` | Sep. 7-Dec. 5 of 2018-2025 | Aligned hourly inflow data in kcfs |
| `soc.csv` | `datetime`, `bon_soc`, `tda_soc`, `jda_soc`, `mcn_soc` | Sep. 7-Dec. 5 of 2018-2025 | Aligned and scaled historical forebay trajectories |

Each processed file should contain 17,280 rows: 90 days x 24 hours x 8 years. Export with `index=False` so the CSV does not contain an extra unnamed index column.

## How `inflow.csv` is generated

1. Load the four raw flow files.
2. Rename the generation, inflow, and spill fields to the common names above.
3. Backward-fill missing values in each raw reservoir dataframe with `bfill()`.
4. Inner-join the four reservoir inflow series on `Date Time`.
5. Retain the 90-day dry-season window from September 7 through December 5 for each complete year from 2018 through 2025.
6. Rename `Date Time` to `datetime` and remove `_kcfs` from the inflow column names.
7. Average duplicate timestamps, construct the complete hourly index for all eight seasons, and reindex to add missing hours.
8. Fill missing hourly values using backward fill and export `inflow.csv`.

The notebook calculates a single minimum and maximum across all four reservoirs and all selected dry-season observations:

$$q_{\min}=\min_{i,t}q_{i,t}, \qquad q_{\max}=\max_{i,t}q_{i,t}.$$

The corresponding normalization would be

$$q^{\mathrm{norm}}_{i,t}=\frac{q_{i,t}-q_{\min}}{q_{\max}-q_{\min}}.$$

However, this transformation is currently commented out in the notebook. Therefore, the present `inflow.csv` contains inflow in kcfs, despite the dataframe name `inflow_data_norm` and the plot label "Normalized inflow."

## How `soc.csv` is generated

1. Load the four raw `*-soc.csv` forebay files and rename each elevation field `soc_ft`.
2. Inner-join all four reservoirs on `Date Time`.
3. Inner-join the result with the timestamps retained in `inflow.csv`, ensuring identical dry-season coverage.
4. Min-max normalize each reservoir's forebay elevation independently over the retained dataset.
5. Scale each normalized series by the corresponding modeled maximum storage:

| Reservoir | Scale $V_i^{\max}$ |
|---|---:|
| Bonneville | 1.5 |
| The Dalles | 1.0 |
| John Day | 1.5 |
| McNary | 4.0 |

For reservoir $i$, the scaled SOC is

$$\mathrm{SOC}_{i,t}=V_i^{\max}\frac{h_{i,t}-\min_t h_{i,t}}{\max_t h_{i,t}-\min_t h_{i,t}},$$

where $h_{i,t}$ is historical forebay elevation in feet. Thus, each output trajectory ranges from 0 to its assigned $V_i^{\max}$ over the retained sample.

6. Remove `_ft` from the column names, average duplicate timestamps, and reindex to the same 17,280-hour index used for inflow.
7. Fill missing hourly values using backward fill and export `soc.csv`.

## Missing-data handling

The notebook uses backward fill, not linear interpolation. A missing observation is replaced by the next available observation. Before filling, the notebook identified 56 missing inflow timestamps and 70 missing SOC timestamps across the 2018-2025 dry-season windows. Re-run these checks whenever the raw USBR files are updated.

## Reproducibility checks

After regenerating the processed files, verify:

```python
assert len(inflow_data_norm) == 17_280
assert len(soc_df_norm) == 17_280
assert not inflow_data_norm.isna().any().any()
assert not soc_df_norm.isna().any().any()
assert inflow_data_norm["datetime"].is_unique
assert soc_df_norm["datetime"].is_unique
assert inflow_data_norm["datetime"].equals(soc_df_norm["datetime"])
```

Recommended exports:

```python
inflow_data_norm.to_csv("inflow.csv", index=False)
soc_df_norm.to_csv("soc.csv", index=False)
```

## Provenance note

The raw files should remain unchanged after download. When the data are refreshed, record the USBR query or download URL, access date, and original date range here so that future versions of `inflow.csv` and `soc.csv` can be reproduced exactly.

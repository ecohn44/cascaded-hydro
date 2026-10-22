"""
fit_inflow_fullyear.py
======================
NEW script. Nothing else is modified -- not the notebook, not the MATLAB.

Purpose
-------
The published fit (`inflow-models 3.ipynb`) is restricted to Sep 7 - Dec 5 by
one line:

    expected = pd.DatetimeIndex([]).append(
        [pd.date_range(f"{y}-09-07", f"{y}-12-05 23:00", freq="h") for y in years])

Everything downstream inherits it. That is the paper's 90-day case-study
window, and it contains zero spill -- but the MATLAB controller runs the full
water year, and its violations are concentrated in May-June. So the DDU
parameters are estimated on the season where the constraint never binds.

This script rebuilds the whole chain over the FULL water year and reports it
next to the Sep-Dec fit.

What the random variable is
---------------------------
Paper eq. (1c):  q_i,k = u_{i-1,k-tau} + s_{i-1,k-tau} + L_i,k
Paper eq. (3) :  q_t   = mu_t + xi_t,  xi ~ N(0, Sigma_t)
Paper eq. (4) :  mu_i,t = a0 + a1 q_{i,t-l} + b1 u_{i-1,t-tau}

These agree once you see that the mean model already carries the routed
release. xi = q - mu is then the part of TOTAL inflow the upstream decision
does not explain -- operationally the incremental-inflow innovation. So
"variance fitted on total inflow" and "variance applied to the local term in
the water balance" are the same object, and no MATLAB change is implied.

That equivalence needs b1 to be identified. In the notebook's 24-lag ridge the
own-lag coefficient is ~1.1 and the release coefficients sum to ~0.02: lagged
inflow already contains last hour's routed release, so persistence absorbs
routing and b1 collapses. This script therefore fits the mean model in the
PAPER'S OWN FORM (eq. 4, three parameters) as the primary specification, and
reports the 24-lag ridge alongside for reference.

2026-09-26 -- SPILL FIX. `{u}_outflow` was Flow-Gen only, with Flow-Spill
read but never added. That is exact inside the Sep7-Dec5 window, where spill
is 0, and wrong over the full water year, where spill runs 40-72% of hours and
54% of the water in May. It made the reconstruction (inflow = outflow + dV/dt)
understate inflow and hid most of the routed release from the upstream
regressors. Outflow is now total release u + s throughout. Effect on the
notebook's own specification (fullyear_ridge24): test R2 0.924 -> 0.982 at JDA,
sigma_DIU 11.00 -> 8.27, and both gamma terms become active. window_ridge24 is
unchanged to two decimals, as it must be.

Outputs (written next to this script):
    fullyear_archx.csv          ARCH-X coefficients, every specification
    fullyear_meanmodel.csv      mean-model coefficients and skill
    fullyear_R_<spec>.csv       residual correlation matrices
    fullyear_report.txt
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
from sklearn.linear_model import RidgeCV, LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

import os
_CANDS = [("/Users/elizacohn/Desktop/cascaded-hydro/simulation-data",
           "/Users/elizacohn/Desktop/Columbia River Data/US ACE Data 2026"),
          ("/mnt/user-data/uploads/cascaded-hydro/simulation-data",
           "/mnt/user-data/uploads/Columbia River Data/US ACE Data 2026")]
SIM, ACE = next((s, a) for s, a in _CANDS if os.path.isdir(s))
OUT = os.path.dirname(os.path.abspath(__file__))

UNITS = ["mcn", "jda", "tda", "bon"]
K = 0.0826446                      # kaf per kcfs-hour
HEADFIT = {"bon": (0.85, 0.08), "tda": (0.73, 0.21),
           "jda": (0.83, 0.12), "mcn": (1.0, 0.0)}
# unit -> (upstream, own lag l, release travel lag tau)  [notebook mid-lags]
REACH = {"jda": ("mcn_outflow", 1, 10),
         "tda": ("jda_outflow", 1, 2),
         "bon": ("tda_outflow", 1, 4)}
RIDGE_LAGS = {"jda": (24, [9, 10, 11]), "tda": (24, [1, 2, 3]), "bon": (24, [3, 4, 5])}
TRAIN = [2018, 2019, 2020, 2021, 2022]
TEST = [2023, 2024, 2025]
CAP = 0.995

LOG = []


def say(s=""):
    print(s)
    LOG.append(str(s))


# --------------------------------------------------------------------------
# 1. rebuild the hourly inflow reconstruction over the full record
# --------------------------------------------------------------------------
def build_recon():
    src = {"bon": "bonneville", "tda": "dalles", "jda": "johnday", "mcn": "mcnary"}
    pre = {"bon": "BON", "tda": "TDA", "jda": "JDA", "mcn": "MCN"}

    fore = None
    for u in UNITS:
        d = pd.read_csv(f"{ACE}/{src[u]}-head.csv")
        d["datetime"] = pd.to_datetime(d["Date Time"], format="%d-%b-%Y %H:%M",
                                       errors="coerce")
        col = f"{pre[u]}.Elev-Forebay.Inst.1Hour.0.CBT-REV [ft]"
        k = d[["datetime", col]].rename(columns={col: f"{u}_forebay_ft"}).dropna()
        fore = k if fore is None else fore.merge(k, on="datetime", how="inner")

    out = None
    for u in UNITS:
        d = pd.read_csv(f"{SIM}/{src[u]}/{src[u]}.csv")
        d["datetime"] = pd.to_datetime(d["Date Time"], errors="coerce")
        g = f"{pre[u]}.Flow-Gen.Ave.1Hour.1Hour.CBT-REV [kcfs]"
        s = f"{pre[u]}.Flow-Spill.Ave.1Hour.1Hour.CBT-REV [kcfs]"
        keep = d[["datetime", g]].rename(columns={g: f"{u}_outflow"})
        # ---- SPILL FIX (2026-09-26) -------------------------------------
        # outflow was Flow-Gen only. Exact in the Sep7-Dec5 window, where
        # spill is 0, but over the full water year spill is 40-72% of hours
        # and up to 54% of the water in May, so gen-only both understates
        # the reconstructed inflow (inflow = outflow + dV/dt) and hides most
        # of the routed release from the upstream regressors.
        # Total release u + s is now used throughout.
        if s in d.columns:
            sp = pd.to_numeric(d[s], errors="coerce").fillna(0.0).values
            keep[f"{u}_spill"] = sp
            keep[f"{u}_outflow"] = keep[f"{u}_outflow"].to_numpy(float) + sp
        # -----------------------------------------------------------------
        out = keep if out is None else out.merge(keep, on="datetime", how="inner")

    df = out.merge(fore, on="datetime", how="inner").sort_values("datetime")
    df = df.dropna(subset=[f"{u}_forebay_ft" for u in UNITS]).reset_index(drop=True)

    rec = pd.DataFrame({"datetime": df["datetime"]})
    for u in UNITS:
        a, b = HEADFIT[u]
        H = df[f"{u}_forebay_ft"]
        Hn = (H - H.min()) / (H.max() - H.min())
        Vn = np.clip((Hn - b) / a, 0, 1)
        Vkaf = Vn.min() + Vn * (Vn.max() - Vn.min())
        rec[f"{u}_inflow"] = df[f"{u}_outflow"] + Vkaf.diff() / K
    for u in UNITS:                                    # bias-correct as the notebook does
        rec[f"{u}_inflow"] -= rec[f"{u}_inflow"].mean() - df[f"{u}_outflow"].mean()
    for c in ["_outflow", "_spill"]:
        for u in UNITS:
            if f"{u}{c}" in df.columns:
                rec[f"{u}{c}"] = df[f"{u}{c}"].values

    rec = rec.dropna().reset_index(drop=True)
    dt = rec["datetime"]
    rec["season_year"] = dt.dt.year - (dt.dt.month <= 4)
    rec["in_window"] = ((dt.dt.month == 9) & (dt.dt.day >= 7)) | \
                       (dt.dt.month.isin([10, 11])) | \
                       ((dt.dt.month == 12) & (dt.dt.day <= 5))
    return rec


# --------------------------------------------------------------------------
# 2. mean models
# --------------------------------------------------------------------------
def mean_paper(d, unit, train_years):
    """Paper eq. (4): mu = a0 + a1 q_{t-l} + b1 u_{up,t-tau}. Three parameters."""
    up, l, tau = REACH[unit]
    t = f"{unit}_inflow"
    g = d.groupby("season_year")
    x = pd.DataFrame({"q_lag": g[t].shift(l), "u_lag": g[up].shift(tau)})
    y = d[t]
    ok = x.notna().all(1) & y.notna()
    tr = ok & d["season_year"].isin(train_years)
    m = LinearRegression().fit(x[tr], y[tr])
    pred = pd.Series(np.nan, index=d.index)
    pred[ok] = np.maximum(0.0, m.predict(x[ok]))
    return pred, dict(a0=float(m.intercept_), a1=float(m.coef_[0]), b1=float(m.coef_[1]))


def mean_ridge(d, unit, train_years):
    """The notebook's specification: p own lags + 3 release lags, RidgeCV."""
    up = REACH[unit][0]
    p, rl = RIDGE_LAGS[unit]
    t = f"{unit}_inflow"
    g = d.groupby("season_year")
    X = pd.DataFrame({f"q{h}": g[t].shift(h) for h in range(1, p + 1)})
    for h in rl:
        X[f"r{h}"] = g[up].shift(h)
    y = d[t]
    ok = X.notna().all(1) & y.notna()
    tr = ok & d["season_year"].isin(train_years)
    m = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-4, 4, 40)))
    m.fit(X[tr], y[tr])
    pred = pd.Series(np.nan, index=d.index)
    pred[ok] = np.maximum(0.0, m.predict(X[ok]))
    sc, rg = m[0], m[-1]
    c = rg.coef_ / sc.scale_
    return pred, dict(a0=float(rg.intercept_ - c @ sc.mean_), a1=float(c[0]),
                      b1=float(c[p:].sum()))


def mcn_mean(d):
    return d.groupby("season_year")["mcn_inflow"].shift(1)


# --------------------------------------------------------------------------
# 3. ARCH-X
# --------------------------------------------------------------------------
def archx_frame(d, unit, resid):
    if unit == "mcn":
        f = pd.DataFrame({"residual": resid, "season_year": d["season_year"],
                          "release": 0.0, "ramp": 0.0})
        return f.dropna()
    up, _, tau = REACH[unit]
    g = d.groupby("season_year")
    f = pd.DataFrame({
        "residual": resid,
        "season_year": d["season_year"],
        "release": g[up].shift(tau),
        "ramp": g[up].diff().groupby(d["season_year"]).shift(tau).abs()})
    return f.dropna()


def fit_archx(tr):
    e = tr["residual"].to_numpy(float)
    v0 = float(np.mean(e ** 2))
    sy = tr["season_year"].to_numpy()
    reset = np.r_[True, sy[1:] != sy[:-1]]
    rng, X = {}, []
    for c in ["release", "ramp"]:
        lo = float(tr[c].min())
        sp = float(tr[c].max() - lo) or 1.0
        X.append((tr[c].to_numpy(float) - lo) / sp)
        rng[c] = (lo, sp)
    X = np.column_stack(X)
    e2 = np.r_[0.0, e[:-1] ** 2]

    def nll(p):
        w, a, gu, gr = p
        s2 = np.where(reset, v0, w + a * e2 + gu * X[:, 0] + gr * X[:, 1])
        if np.any(~np.isfinite(s2)) or np.any(s2 <= 0):
            return 1e20
        return -(-0.5 * (np.log(2 * np.pi) + np.log(s2) + e ** 2 / s2))[~reset].sum()

    best = min((minimize(nll, s, method="L-BFGS-B",
                         bounds=[(1e-8, None), (0, CAP), (0, None), (0, None)],
                         options={"maxiter": 5000, "ftol": 1e-8})
                for s in ([0.8 * v0, 0.10, 0.05 * v0, 0.05 * v0],
                          [0.5 * v0, 0.25, 0.10 * v0, 0.10 * v0])),
               key=lambda f: f.fun)
    w, a, gu, gr = best.x
    return dict(omega=w, alpha=a, beta=0.0, gamma_release=gu, gamma_ramp=gr,
                x_ranges=rng, v0=v0, converged=bool(best.success))


def apply_archx(te, m):
    e = te["residual"].to_numpy(float)
    z = []
    for c in ["release", "ramp"]:
        lo, sp = m["x_ranges"][c]
        z.append(np.maximum((te[c].to_numpy(float) - lo) / sp, 0.0))
    sy = te["season_year"].to_numpy()
    reset = np.r_[True, sy[1:] != sy[:-1]]
    e2 = np.r_[0.0, e[:-1] ** 2]
    v = np.where(reset, m["v0"],
                 m["omega"] + m["alpha"] * e2 + m["gamma_release"] * z[0]
                 + m["gamma_ramp"] * z[1])
    v = np.maximum(v, 1e-12)
    return np.sqrt(v), e / np.sqrt(v), reset


# --------------------------------------------------------------------------
def run(rec, spec, mean_fn, window_only):
    d = rec[rec["in_window"]].reset_index(drop=True) if window_only else rec
    rows, mrows, resid = [], [], {}

    for unit in UNITS:
        if unit == "mcn":
            pred = mcn_mean(d)
            coef = dict(a0=0.0, a1=1.0, b1=0.0)
        else:
            pred, coef = mean_fn(d, unit, TRAIN)
        r = d[f"{unit}_inflow"] - pred
        resid[unit] = r
        msk = r.notna()
        te = msk & d["season_year"].isin(TEST)
        tr = msk & d["season_year"].isin(TRAIN)
        mrows.append(dict(spec=spec, unit=unit, **coef,
                          sigma_DIU=float(np.sqrt(np.mean(r[tr] ** 2))),
                          test_R2=float(r2_score(d[f"{unit}_inflow"][te], pred[te]))
                          if te.sum() > 10 else np.nan,
                          n=int(msk.sum())))

        F = archx_frame(d, unit, r)
        Ftr = F[F.season_year.isin(TRAIN)]
        Fte = F[F.season_year.isin(TEST)]
        if len(Ftr) < 500 or len(Fte) < 100:
            continue
        m = fit_archx(Ftr)
        if unit == "mcn":
            m["gamma_release"] = m["gamma_ramp"] = 0.0
        sg, zz, rs = apply_archx(Fte, m)
        rows.append(dict(spec=spec, unit=unit, omega=m["omega"], alpha=m["alpha"],
                         beta=0.0, gamma_release=m["gamma_release"],
                         gamma_ramp=m["gamma_ramp"],
                         sigma_DIU=mrows[-1]["sigma_DIU"],
                         release_min=m["x_ranges"]["release"][0],
                         release_range=m["x_ranges"]["release"][1],
                         ramp_min=m["x_ranges"]["ramp"][0],
                         ramp_range=m["x_ranges"]["ramp"][1],
                         n_train=len(Ftr), n_test=len(Fte),
                         sigma_mean=float(sg[~rs].mean()),
                         sigma_p5=float(np.percentile(sg[~rs], 5)),
                         sigma_p95=float(np.percentile(sg[~rs], 95)),
                         coverage_95=float((np.abs(zz[~rs]) <= 1.959964).mean()),
                         converged=m["converged"]))

    E = pd.DataFrame(resid)
    E = E[d["season_year"].isin(TEST).values].dropna()
    return pd.DataFrame(rows), pd.DataFrame(mrows), E.corr().loc[UNITS, UNITS]


# --------------------------------------------------------------------------
# 4. extra exports for the Fig. 2 redesign (v2)
# --------------------------------------------------------------------------
def export_panel_data(rec):
    """Panel (b): predicted vs. observed on held-out years, all 4 dams, under
    the fullyear_ridge24 spec already hardcoded in make_fig2.m -- so the
    'four heavy-tailed distributions' overlay is the exact same model as the
    scatter it replaces, just extended from TDA to all four units.

    Panels (c)/(d): the full per-timestep TDA series (residual, upstream
    release/ramp, and the fitted ARCH-X sigma path) so the plotting code can
    build violin bins and a time-series band plot without refitting anything.
    Nothing here changes what run() already computed or wrote; this only
    dumps intermediate arrays that main()'s loop discards after aggregating.
    """
    d = rec  # fullyear_ridge24 is the unwindowed (window_only=False) spec
    b_rows = []
    tda_resid = None
    v0_rows = []
    for unit in UNITS:
        if unit == "mcn":
            pred = mcn_mean(d)
        else:
            pred, _ = mean_ridge(d, unit, TRAIN)
        r = d[f"{unit}_inflow"] - pred
        te = r.notna() & d["season_year"].isin(TEST)
        b_rows.append(pd.DataFrame({
            "unit": unit, "datetime": d["datetime"][te],
            "observed": d[f"{unit}_inflow"][te], "predicted": pred[te]}))
        if unit == "tda":
            tda_resid = r
        Fu = archx_frame(d, unit, r)
        Fu_tr = Fu[Fu.season_year.isin(TRAIN)]
        v0_rows.append(dict(unit=unit, v0=float(np.mean(Fu_tr["residual"].to_numpy(float) ** 2))))

    Bdf = pd.concat(b_rows, ignore_index=True)
    Bdf.to_csv(f"{OUT}/fullyear_panelb_test.csv", index=False)
    pd.DataFrame(v0_rows).to_csv(f"{OUT}/fullyear_v0.csv", index=False)

    F = archx_frame(d, "tda", tda_resid)
    F["datetime"] = d.loc[F.index, "datetime"].values
    F["in_window"] = d.loc[F.index, "in_window"].values
    F["is_train"] = F["season_year"].isin(TRAIN)
    F["is_test"] = F["season_year"].isin(TEST)
    m = fit_archx(F[F.is_train])
    sg, zz, reset = apply_archx(F, m)
    F["sigma_archx"] = sg
    F["z"] = zz
    F["season_reset"] = reset
    F["sigma_DIU"] = float(np.sqrt(np.mean(tda_resid[d["season_year"].isin(TRAIN)] ** 2)))
    F.to_csv(f"{OUT}/fullyear_paneld_tda_series.csv", index=False)

    # full per-unit ridge coefficients (fullyear_ridge24 spec) so the MATLAB
    # port can hardcode all 4 dams the same way TDA already is -- mean_ridge()
    # only returns the compressed a0/a1/b1 for the report, this recovers the
    # raw p own-lag + release-lag coefficient vector it actually fits with.
    coef_rows = []
    for unit in UNITS:
        if unit == "mcn":
            coef_rows.append(dict(unit=unit, intercept=0.0, own_lags=list(range(1, 2)),
                                   own_coefs=[1.0], release_lags=[], release_coefs=[]))
            continue
        up = REACH[unit][0]
        p, rl = RIDGE_LAGS[unit]
        t = f"{unit}_inflow"
        g = d.groupby("season_year")
        X = pd.DataFrame({f"q{h}": g[t].shift(h) for h in range(1, p + 1)})
        for h in rl:
            X[f"r{h}"] = g[up].shift(h)
        y = d[t]
        ok = X.notna().all(1) & y.notna()
        tr = ok & d["season_year"].isin(TRAIN)
        mm = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-4, 4, 40)))
        mm.fit(X[tr], y[tr])
        sc, rg = mm[0], mm[-1]
        c = rg.coef_ / sc.scale_
        b0 = float(rg.intercept_ - c @ sc.mean_)
        coef_rows.append(dict(unit=unit, intercept=b0, own_lags=list(range(1, p + 1)),
                               own_coefs=c[:p].tolist(), release_lags=rl,
                               release_coefs=c[p:].tolist()))
    import json
    with open(f"{OUT}/fullyear_ridge_coeffs_full.json", "w") as fh:
        json.dump(coef_rows, fh, indent=1)

    # panel (a): the full reconstructed hourly series for all 4 units, straight
    # from build_recon() -- so panel (a) can be redrawn unchanged, or checked
    # against whichever unit it was actually plotting.
    Hdf = d[["datetime", "season_year"] + [f"{u}_inflow" for u in UNITS]].copy()
    Hdf.to_csv(f"{OUT}/fullyear_hydrograph_all_units.csv", index=False)

    # panel (c) diagnostic: how well does the fitted ARCH-X *Gaussian* actually
    # capture the conditional skew, bin by bin? Two measures per bin:
    #  (1) Jarque-Bera test on the standardized residual z = e/sigma_archx
    #      (jointly tests skew and excess kurtosis against the Gaussian null)
    #  (2) one-sided empirical vs. nominal tail exceedance at each eps the
    #      paper actually uses (0.01/0.05/0.10/0.20) -- this is the number
    #      that matters operationally: does z_{1-eps}*sigma_archx actually
    #      bound (1-eps) of the mass on each side, or does skew make one
    #      tail worse than the nominal eps promises?
    from scipy.stats import jarque_bera, skew as _skew, kurtosis as _kurt, norm as _norm
    rel = F.release.to_numpy(float)
    z = F.z.to_numpy(float)
    qs = np.quantile(rel, np.linspace(0, 1, 8))
    qs[0] -= 1e-6
    idx = np.clip(np.digitize(rel, qs) - 1, 0, len(qs) - 2)
    EPS = [0.01, 0.05, 0.10, 0.20]
    diag_rows = []
    for b in range(len(qs) - 1):
        msk = idx == b
        if msk.sum() < 200:
            continue
        zb = z[msk]
        jb_stat, jb_p = jarque_bera(zb)
        row = dict(bin=b, release_center=float(rel[msk].mean()), n=int(msk.sum()),
                   skew=float(_skew(zb)), excess_kurtosis=float(_kurt(zb)),
                   jarque_bera_stat=float(jb_stat), jarque_bera_p=float(jb_p))
        for e in EPS:
            zc = _norm.ppf(1 - e)
            row[f"upper_exceed_nominal_{e}"] = e
            row[f"upper_exceed_empirical_{e}"] = float((zb > zc).mean())
            row[f"lower_exceed_nominal_{e}"] = e
            row[f"lower_exceed_empirical_{e}"] = float((zb < -zc).mean())
        diag_rows.append(row)
    Cdiag = pd.DataFrame(diag_rows)
    Cdiag.to_csv(f"{OUT}/fullyear_panelc_normality_diagnostics.csv", index=False)

    say(f"wrote fullyear_panelb_test.csv ({len(Bdf)} rows, 4 units, held-out years), "
        f"fullyear_paneld_tda_series.csv ({len(F)} rows, TDA full record), "
        f"fullyear_hydrograph_all_units.csv ({len(Hdf)} rows), "
        f"fullyear_panelc_normality_diagnostics.csv ({len(Cdiag)} bins), "
        f"fullyear_ridge_coeffs_full.json ({len(coef_rows)} units)")


def main():
    rec = build_recon()
    say(f"reconstruction rebuilt: {len(rec)} hourly rows, "
        f"{rec.datetime.min()} -> {rec.datetime.max()}")
    say(f"  months present: {sorted(rec.datetime.dt.month.unique().tolist())}")
    say(f"  Sep7-Dec5 window is {100*rec['in_window'].mean():.1f}% of it")
    if "tda_spill" in rec:
        say(f"  TDA spill>0: {100*(rec.tda_spill>0).mean():.1f}% full year vs "
            f"{100*(rec.loc[rec.in_window,'tda_spill']>0).mean():.1f}% in the window")
    say()

    specs = [("window_paperform", mean_paper, True),
             ("fullyear_paperform", mean_paper, False),
             ("window_ridge24", mean_ridge, True),
             ("fullyear_ridge24", mean_ridge, False)]

    A, M, Rs = [], [], {}
    for name, fn, win in specs:
        say(f"===== {name}")
        a, m, R = run(rec, name, fn, win)
        A.append(a); M.append(m); Rs[name] = R
        if len(a):
            say(a[["unit", "omega", "alpha", "gamma_release", "gamma_ramp",
                   "sigma_DIU", "sigma_mean", "coverage_95", "n_train"]]
                .round(3).to_string(index=False))
        say(m[["unit", "a0", "a1", "b1", "sigma_DIU", "test_R2"]].round(4)
            .to_string(index=False))
        say("  R (test residuals):")
        say("  " + R.round(4).to_string().replace("\n", "\n  "))
        say()

    A = pd.concat(A, ignore_index=True); M = pd.concat(M, ignore_index=True)
    A.to_csv(f"{OUT}/fullyear_archx.csv", index=False)
    M.to_csv(f"{OUT}/fullyear_meanmodel.csv", index=False)
    for k, R in Rs.items():
        R.to_csv(f"{OUT}/fullyear_R_{k}.csv")

    say("===== what the Sep-Dec window was costing (paper-form mean model)")
    w = A[A.spec == "window_paperform"].set_index("unit")
    f = A[A.spec == "fullyear_paperform"].set_index("unit")
    for u in UNITS:
        if u in w.index and u in f.index:
            say(f"  {u.upper():4s}  sigma_DIU {w.loc[u,'sigma_DIU']:6.2f} -> {f.loc[u,'sigma_DIU']:6.2f}"
                f"   g_rel {w.loc[u,'gamma_release']:8.2f} -> {f.loc[u,'gamma_release']:8.2f}"
                f"   g_ramp {w.loc[u,'gamma_ramp']:8.2f} -> {f.loc[u,'gamma_ramp']:8.2f}"
                f"   cover95 {w.loc[u,'coverage_95']:.3f} -> {f.loc[u,'coverage_95']:.3f}")

    export_panel_data(rec)

    with open(f"{OUT}/fullyear_report.txt", "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    say()
    say(f"wrote fullyear_archx.csv, fullyear_meanmodel.csv, fullyear_R_*.csv, "
        f"fullyear_panelb_test.csv, fullyear_paneld_tda_series.csv to {OUT}")


if __name__ == "__main__":
    main()

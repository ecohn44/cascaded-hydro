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

SIM = "/mnt/user-data/uploads/cascaded-hydro/simulation-data"
ACE = "/mnt/user-data/uploads/Columbia River Data/US ACE Data 2026"
OUT = "/home/claude/fit"

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
        if s in d.columns:
            keep[f"{u}_spill"] = pd.to_numeric(d[s], errors="coerce").values
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

    with open(f"{OUT}/fullyear_report.txt", "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    say()
    say(f"wrote fullyear_archx.csv, fullyear_meanmodel.csv, fullyear_R_*.csv to {OUT}")


if __name__ == "__main__":
    main()

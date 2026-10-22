"""Raw variance diagnostics for the ARCH-X specification (Fig. 2 family).

Replicates make_fig2.m's reconstruction and mean model EXACTLY -- same
forebay->volume shape, same ridge24 coefficients, same season-year reset --
then stops before the ARCH-X sigma path and looks at the residuals directly.

Everything here is TDA (the dam Fig. 2c/2e use), test years WY2023-25.
"""

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures_out'
OUT.mkdir(exist_ok=True)

import numpy as np, pandas as pd

ACE = str(ROOT) + '/data/raw/usace-forebay/'
SIM = str(ROOT) + '/data/raw/usace-flow/'
units = ['mcn', 'jda', 'tda', 'bon']
proj = ['mcnary', 'johnday', 'dalles', 'bonneville']
pre = ['MCN', 'JDA', 'TDA', 'BON']
headA = [1.00, 0.83, 0.73, 0.85]
headB = [0.00, 0.12, 0.21, 0.08]
K = 0.0826446
CMS = 28.3168
TRAIN, TEST = range(2018, 2023), range(2023, 2026)

MDL = {
 'jda': dict(up='mcn', tau=10, B0=0.24225651729057063, QLAG=list(range(1, 25)),
   QCOEF=[1.12074676,-0.23143458,-0.00912926,-0.00469043,-0.00495693,-0.02671501,
          -0.02952234,-0.02662365,0.06508141,0.00061729,0.01556349,0.01556629,
          -0.01286978,-0.01749457,-0.01029470,-0.05393533,0.02530521,0.01839794,
          0.02234139,0.00845066,0.01049564,0.03165186,0.13894447,-0.07860834],
   RLAG=[9,10,11], RCOEF=[0.01957648,0.00695373,0.00548536],
   omega=39.86571805966577, alpha=0.2252424666315235, g_rel=51.787117733902754,
   g_ramp=23.634242128892804, v0=68.3960856481441,
   rel_min=47.9, rel_rng=469.5, rmp_min=0.0, rmp_rng=74.39999999999998,
   sigma_DIU=8.270192600426189),
 'tda': dict(up='jda', tau=2, B0=0.5197505746442914, QLAG=list(range(1, 25)),
   QCOEF=[0.89714400,-0.15788527,0.02877379,-0.01188400,-0.01417420,-0.03312762,
          -0.01455826,-0.00084958,0.05319990,0.00790728,0.01585623,-0.00911670,
          0.00331579,-0.01981391,-0.00272361,-0.05251301,0.02447719,0.01244897,
          0.01294586,0.01514867,0.00214110,0.02506731,0.13084201,-0.08391526],
   RLAG=[1,2,3], RCOEF=[0.36615223,-0.21533754,0.00757301],
   omega=35.93363590660807, alpha=0.1736814331710667, g_rel=61.12944733824949,
   g_ramp=14.634112398619068, v0=60.6423311623629,
   rel_min=50.8, rel_rng=466.3, rmp_min=0.0, rmp_rng=101.09999999999997,
   sigma_DIU=7.787318611843418),
 'bon': dict(up='tda', tau=4, B0=0.3063089512114914, QLAG=list(range(1, 25)),
   QCOEF=[1.23337297,-0.30097188,0.00436717,-0.00683676,-0.01495385,-0.02616897,
          0.00246222,-0.01315426,0.03047101,0.01889004,0.00052087,0.00752304,
          -0.00382333,-0.00424788,-0.02447322,-0.00501978,0.01339631,0.00644078,
          0.01249213,0.00458486,0.01846936,0.02785399,0.09847428,-0.10171162],
   RLAG=[3,4,5], RCOEF=[0.02864034,-0.00915402,0.00184186],
   omega=10.521644577201169, alpha=0.2451478751985036, g_rel=34.47408216689823,
   g_ramp=0.0, v0=24.914305665223797,
   rel_min=23.4, rel_rng=468.6, rmp_min=0.0, rmp_rng=118.1,
   sigma_DIU=4.991423210390379),
}


def load():
    fb, fl = None, None
    for i in range(4):
        h = pd.read_csv(ACE + proj[i] + '-head.csv')
        t = pd.to_datetime(h['Date Time'], format='%d-%b-%Y %H:%M', errors='coerce')
        col = '%s.Elev-Forebay.Inst.1Hour.0.CBT-REV [ft]' % pre[i]
        k = pd.DataFrame({'dt': t, units[i] + '_fb': pd.to_numeric(h[col], errors='coerce')}).dropna()
        fb = k if fb is None else fb.merge(k, on='dt', how='inner')

        d = pd.read_csv(SIM + proj[i] + '/' + proj[i] + '.csv')
        t = pd.to_datetime(d['Date Time'], errors='coerce')
        gen = pd.to_numeric(d['%s.Flow-Gen.Ave.1Hour.1Hour.CBT-REV [kcfs]' % pre[i]], errors='coerce')
        spl = pd.to_numeric(d['%s.Flow-Spill.Ave.1Hour.1Hour.CBT-REV [kcfs]' % pre[i]], errors='coerce').fillna(0)
        k = pd.DataFrame({'dt': t, units[i] + '_out': gen + spl}).dropna(subset=['dt'])
        fl = k if fl is None else fl.merge(k, on='dt', how='inner')
    M = fl.merge(fb, on='dt', how='inner').sort_values('dt').reset_index(drop=True)
    M = M.dropna(subset=[u + '_fb' for u in units])
    for i in range(4):
        H = M[units[i] + '_fb'].values
        Hn = (H - H.min()) / (H.max() - H.min())
        Vn = np.clip((Hn - headB[i]) / headA[i], 0, 1)
        Vk = Vn.min() + Vn * (Vn.max() - Vn.min())
        q = M[units[i] + '_out'].values + np.r_[np.nan, np.diff(Vk)] / K
        q = q - (np.nanmean(q) - np.nanmean(M[units[i] + '_out'].values))
        M[units[i] + '_q'] = q
    M['sy'] = M.dt.dt.year - (M.dt.dt.month <= 4).astype(int)
    return M.reset_index(drop=True)


def shift_season(x, h, new):
    y = np.r_[np.full(h, np.nan), x[:-h]] if h > 0 else x.copy()
    bad = np.zeros(len(x), bool)
    idx = np.flatnonzero(new)
    for s in idx:
        bad[s:s + h] = True
    y[bad] = np.nan
    return y


def residuals(M, u):
    m = MDL[u]
    q = M[u + '_q'].values
    new = np.r_[True, np.diff(M.sy.values) != 0]
    X = np.column_stack(
        [shift_season(q, h, new) for h in m['QLAG']] +
        [shift_season(M[m['up'] + '_out'].values, j, new) for j in m['RLAG']])
    pred = np.maximum(0, m['B0'] + X @ np.r_[m['QCOEF'], m['RCOEF']])
    e = q - pred
    uup = M[m['up'] + '_out'].values
    rel = shift_season(uup, m['tau'], new)
    ramp = np.abs(shift_season(np.r_[np.nan, np.diff(uup)], m['tau'], new))
    gd = ~np.isnan(e) & ~np.isnan(rel) & ~np.isnan(ramp)
    E, R, P, SY = e[gd], rel[gd], ramp[gd], M.sy.values[gd]
    relN = np.maximum((R - m['rel_min']) / m['rel_rng'], 0)
    rmpN = np.maximum((P - m['rmp_min']) / m['rmp_rng'], 0)
    rst = np.r_[True, np.diff(SY) != 0]
    e2 = np.r_[0, E[:-1] ** 2]
    s2 = m['omega'] + m['alpha'] * e2 + m['g_rel'] * relN + m['g_ramp'] * rmpN
    s2[rst] = m['v0']
    return dict(E=E, rel=R, ramp=P, SY=SY, sig=np.sqrt(np.maximum(s2, 1e-12)))


if __name__ == '__main__':
    M = load()
    print('rows', len(M), M.dt.min(), '->', M.dt.max())
    for u in ('jda', 'tda', 'bon'):
        r = residuals(M, u)
        isT = np.isin(r['SY'], list(TEST))
        E = r['E'][isT] * CMS
        print('%s  n_test %5d   resid sd %6.1f   ARCH-X sigma mean %6.1f   DIU %6.1f  (m3/s)'
              % (u, isT.sum(), E.std(ddof=1), r['sig'][isT].mean() * CMS,
                 MDL[u]['sigma_DIU'] * CMS))

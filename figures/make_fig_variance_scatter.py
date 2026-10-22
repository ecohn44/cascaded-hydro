"""Raw forecast-error cloud against upstream release, and what a variance
model can actually extract from it.  The Dalles, test years WY2023-25.

Left   every hourly residual, one dot, against effective upstream release
       (the same rel + C_ramp*ramp combination Fig. 2c/2e use). The binned
       EMPIRICAL +/-1 sd envelope is drawn over it, with the constant DIU
       band and the ARCH-X mean band for comparison.
Right  the same binned empirical sd on its own, against the two model
       sigmas. This is the honest version of Fig. 2e.

WHY THIS PANEL EXISTS. Fig. 2e as written bins St = TDA.sig -- the FITTED
ARCH-X sigma -- calls those points "Observed", and then draws a straight
least-squares line through them and calls that line "DDU". Both series are
the model; the empirical spread of the residuals never enters. What is
plotted here instead is std(residual) inside each release bin, which is the
quantity the model is supposed to reproduce.

Residuals, the mean model and the ARCH-X path are reproduced from
make_fig2.m exactly (see variance_probe.py); nothing is refitted.

Output: figures_out/fig_variance_scatter.{pdf,png}
"""

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures_out'
OUT.mkdir(exist_ok=True)

import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from variance_probe import load, residuals, MDL, TEST, CMS

NBIN = 12
M = load()
r = residuals(M, 'tda')
m = MDL['tda']

isT = np.isin(r['SY'], list(TEST))
C = (m['g_ramp'] * m['rel_rng']) / (m['g_rel'] * m['rmp_rng'])
U = (r['rel'] + C * (r['ramp'] - m['rmp_min']))[isT] * CMS / 1000
E = r['E'][isT] * CMS
S = r['sig'][isT] * CMS
DIU = m['sigma_DIU'] * CMS

ed = np.unique(np.quantile(U, np.linspace(0, 1, NBIN + 1)))
ed[0] -= 1e-6
ctr, emp, arx = [], [], []
for b in range(len(ed) - 1):
    s = (U >= ed[b]) & (U < ed[b + 1])
    if s.sum() < 100:
        continue
    ctr.append(U[s].mean())
    emp.append(E[s].std(ddof=1))
    arx.append(S[s].mean())
ctr, emp, arx = map(np.array, (ctr, emp, arx))

DIU_C = np.array([0, 114, 178]) / 255
DDU_C = np.array([213, 94, 0]) / 255
INK = np.array([26, 26, 26]) / 255
RULE = '#dcdcda'
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 8,
    'axes.labelsize': 8, 'legend.fontsize': 7.0,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'axes.edgecolor': '#666', 'axes.linewidth': .6,
    'xtick.color': '#666', 'ytick.color': '#666', 'axes.labelcolor': 'black',
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white'})

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.16, 2.9))

# ---- left: the cloud
a1.axhline(0, color='#bbb', lw=.5, zorder=1)
a1.plot(U, E, '.', ms=1.0, color=INK, alpha=.10, mew=0, zorder=2,
        rasterized=True, label='Hourly forecast error')
a1.plot(ctr, emp, '-', color=DDU_C, lw=1.4, zorder=5,
        label=r'Empirical $\pm1\,\sigma$')
a1.plot(ctr, -emp, '-', color=DDU_C, lw=1.4, zorder=5)
a1.plot([U.min(), U.max()], [DIU] * 2, '--', color=DIU_C, lw=1.3, zorder=4,
        label=r'DIU $\pm1\,\sigma$ (constant)')
a1.plot([U.min(), U.max()], [-DIU] * 2, '--', color=DIU_C, lw=1.3, zorder=4)
a1.set_xlabel(r'Effective upstream release ($10^3$ m$^3$/s)')
a1.set_ylabel(r'Forecast error (m$^3$/s)')
a1.set_ylim(-1100, 1100)
a1.legend(frameon=False, loc='upper left', handletextpad=.5,
          borderaxespad=.2, labelspacing=.3, markerscale=4)

# ---- right: binned sd, the corrected Fig. 2e
a2.plot(ctr, emp, 'o', ms=4.2, mfc=INK, mec=INK, zorder=5,
        label='Empirical (binned)')
a2.plot(ctr, arx, '-', color=DDU_C, lw=1.8, zorder=4, label='ARCH-X')
a2.plot([ctr.min(), ctr.max()], [DIU] * 2, '--', color=DIU_C, lw=1.8,
        zorder=3, label='DIU')
a2.set_xlabel(r'Effective upstream release ($10^3$ m$^3$/s)')
a2.set_ylabel(r'Std. dev. of forecast error (m$^3$/s)')
a2.legend(frameon=False, loc='upper left', handletextpad=.5,
          borderaxespad=.2, labelspacing=.3)

for ax in (a1, a2):
    ax.grid(color=RULE, lw=.5)
    ax.set_axisbelow(True)

fig.tight_layout(pad=.3, w_pad=1.6)
for ext in ('pdf', 'png'):
    fig.savefig(str(OUT) + '/fig_variance_scatter.' + ext, dpi=600,
                bbox_inches='tight', facecolor='white')

print('bin    n      U_ctr   empirical sd   ARCH-X   DIU')
for b in range(len(ed) - 1):
    s = (U >= ed[b]) & (U < ed[b + 1])
    if s.sum() < 100:
        continue
    print('%3d %6d   %7.2f   %10.1f  %8.1f  %6.1f'
          % (b + 1, s.sum(), U[s].mean(), E[s].std(ddof=1), S[s].mean(), DIU))

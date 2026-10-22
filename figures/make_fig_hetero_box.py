"""Forecast error is heteroskedastic -- the submitted box panel, with the
dispersion statistics the chance constraint actually uses drawn over it.

The Dalles, test years WY2023-25, 11 equal-count bins of effective upstream
release (the same rel + C_ramp*ramp combination as Fig. 2c/2e).

BOXES are the distribution of |forecast error| in each bin: median, IQR, and
whiskers at p10/p90 rather than the Tukey 1.5*IQR rule, so the tail is drawn
at a stated quantile instead of wherever the fences happen to land.

LINES are standard deviations: the empirical sd of the (signed) error in each
bin, the ARCH-X sigma the controller actually uses, and the constant DIU sigma.

WHY BOTH. The box centres and the tails disagree, and the disagreement is the
point. Median |error| rises 53 -> 132 m^3/s over the first six bins and then
stops (132, 130, 132, 119, 124, 125). p90 does not stop: 236 -> 503. Neither
does the sd: 152 -> 294. At high release the middle of the error distribution
stops widening while the tail keeps going, i.e. the distribution gets heavier,
not just broader. A panel showing medians alone reads as saturation and
understates the mechanism, because the chance constraint keys on sigma.

Residuals, the mean model and the ARCH-X path are reproduced from make_fig2.m
exactly (see variance_probe.py); nothing is refitted here.

Output: figures_out/fig_hetero_box.{pdf,png}
"""

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures_out'
OUT.mkdir(exist_ok=True)

import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from variance_probe import load, residuals, MDL, TEST, CMS

NBIN = 11
M = load()
r = residuals(M, 'tda')
m = MDL['tda']

isT = np.isin(r['SY'], list(TEST))
C = (m['g_ramp'] * m['rel_rng']) / (m['g_rel'] * m['rmp_rng'])
U = (r['rel'] + C * (r['ramp'] - m['rmp_min']))[isT] * CMS / 1000
E = r['E'][isT] * CMS
A = np.abs(E)
S = r['sig'][isT] * CMS
DIU = m['sigma_DIU'] * CMS

ed = np.unique(np.quantile(U, np.linspace(0, 1, NBIN + 1)))
ed[0] -= 1e-6
groups, ctr, sd, arx = [], [], [], []
for b in range(len(ed) - 1):
    s = (U >= ed[b]) & (U < ed[b + 1])
    if s.sum() < 100:
        continue
    groups.append(A[s])
    ctr.append(U[s].mean())
    sd.append(E[s].std(ddof=1))
    arx.append(S[s].mean())
ctr, sd, arx = map(np.array, (ctr, sd, arx))
x = np.arange(1, len(groups) + 1)

BOX = '#b9d3e6'          # submitted panel's light blue
MED = '#c1502e'          # submitted panel's orange median rule
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

fig, ax = plt.subplots(figsize=(3.5, 3.0))

bp = ax.boxplot(groups, positions=x, widths=.62, whis=(10, 90), showfliers=False,
                patch_artist=True, zorder=3)
for p in bp['boxes']:
    p.set(facecolor=BOX, edgecolor='#5b7c91', linewidth=.7)
for p in bp['medians']:
    p.set(color=MED, linewidth=1.3)
for k in ('whiskers', 'caps'):
    for p in bp[k]:
        p.set(color='#5b7c91', linewidth=.7)

ax.plot(x, sd, '-', color=DDU_C, lw=1.5, zorder=5,
        label='Empirical std. dev.')
ax.plot(x, arx, '--', color=INK, lw=1.2, zorder=5, label=r'ARCH-X $\sigma$')
ax.axhline(DIU, ls='--', color=DIU_C, lw=1.3, zorder=4, label=r'DIU $\sigma$')

ax.set_xticks(x)
ax.set_xticklabels(['%.1f' % c for c in ctr])
ax.set_xlabel(r'Effective upstream release $\tilde u_{i-1,t-\tau_i}$'
              r'  ($10^3$ m$^3$/s)')
ax.set_ylabel(r'$|$Forecast error$|$  (m$^3$/s)')
ax.set_ylim(0, 560)
ax.grid(axis='y', color=RULE, lw=.5)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc='upper left', handletextpad=.5,
          borderaxespad=.2, labelspacing=.3)

fig.tight_layout(pad=.3)
for ext in ('pdf', 'png'):
    fig.savefig(str(OUT) + '/fig_hetero_box.' + ext, dpi=600,
                bbox_inches='tight', facecolor='white')

print(' U_ctr      n   med|e|    p75    p90     sd   ARCH-X')
for g, c, s_, a in zip(groups, ctr, sd, arx):
    print('%6.1f %6d %7.0f %6.0f %6.0f %6.0f %8.0f'
          % (c, len(g), np.median(g), np.percentile(g, 75),
             np.percentile(g, 90), s_, a))

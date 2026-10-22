"""Selecting theta from the generation / violation trade-off.

One point per theta: annual generation against integrated violation, DDU, P20,
eps = 0.05, Bonferroni. Marker is the median across scenarios; the shaded
envelope is the p40-p60 band, interpolated along the trajectory (PCHIP, so it
stays monotone between knots and does not overshoot).

The good corner is bottom RIGHT: more energy, less violation.

SAMPLE SIZE IS NOT UNIFORM. theta = 200 to 1400 were run on six scenarios;
100, 1600, 3200 and 5000 remain at three. Filled markers are n = 6, open
markers n = 3. Band widths are therefore not comparable between the two sets,
for reasons that have nothing to do with theta -- say so in the caption.

theta = 0, 10, 25 and 50 are omitted; they sit on top of one another at the
untracked end and only crowd the panel.
"""

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures_out'
OUT.mkdir(exist_ok=True)

import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

B = str(ROOT) + '/results/runs/'
TH = [100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 3200, 5000]
PICK = 400
NMAX = 10       # use every scenario a theta has

g, v, vlo, vhi, nn = [], [], [], [], []
for t in TH:
    d = pd.read_csv(B + 'BPA_TH7_P20_TH%07d/summary_BPA_TH7_P20_TH%07d.csv' % (t, t))
    d = d[(d.policy == 'DDU') & (d.scen <= NMAX)]   # keep the panel BALANCED
    e, i = d.energy_MWh.values / 1e6, d.lower_ivi_kaf.values
    g.append(np.median(e)); v.append(np.median(i)); nn.append(len(e))
    vlo.append(np.percentile(i, 40)); vhi.append(np.percentile(i, 60))
g, v, vlo, vhi, nn = map(np.array, (g, v, vlo, vhi, nn))

INK2, RULE = '#555', '#dcdcda'
C, CH = '#3f7d2a', '#b5342a'
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 8,
    'axes.labelsize': 8, 'legend.fontsize': 7.0,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'axes.edgecolor': '#666', 'axes.linewidth': .6,
    'xtick.color': '#666', 'ytick.color': '#666', 'axes.labelcolor': 'black',
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white'})

fig, ax = plt.subplots(figsize=(3.5, 2.9))

o = np.argsort(g)
xs = np.linspace(g[o][0], g[o][-1], 400)
ax.fill_between(xs, PchipInterpolator(g[o], vlo[o])(xs),
                PchipInterpolator(g[o], vhi[o])(xs),
                color=C, alpha=.15, lw=0, zorder=1)
ax.plot(xs, PchipInterpolator(g[o], v[o])(xs), '-', color=C, lw=1.0,
        alpha=.75, zorder=3)

# One entry per sample size PRESENT, so the legend never carries a dead row.
for k, (ms, fc, mew) in {10: (5.0, C, 1.0), 6: (4.2, C, 1.0),
                         3: (4.2, 'white', 1.1)}.items():
    m = nn == k
    if m.any():
        ax.plot(g[m], v[m], 'o', ms=ms, mfc=fc, mec=C, mew=mew, zorder=4,
                label='$n=%d$ scenarios' % k)

k = TH.index(PICK)
ax.plot(g[k], v[k], 'o', ms=6.4, mfc=CH, mec=CH, zorder=5)

OFF = {100: (4, 6), 200: (5, 5), 400: (7, -12), 600: (5, 5), 800: (-9, -13),
       1000: (1, -13), 1200: (-2, 8), 1400: (5, 5), 1600: (6, 4),
       3200: (5, -13), 5000: (3, 9)}
for i, t in enumerate(TH):
    ax.annotate(r'$\theta=%d$' % t if t in (100, 5000) else r'$%d$' % t,
                (g[i], v[i]), textcoords='offset points', xytext=OFF[t],
                fontsize=7.0, color=CH if t == PICK else INK2)

ax.margins(y=.12)
ax.set_xlabel('Annual generation (TWh)')
ax.set_ylabel('Integrated violation index (kaf)')
ax.grid(color=RULE, lw=.5)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc='upper left', handletextpad=.4, borderaxespad=.2)

fig.tight_layout(pad=.25)
for ext in ('pdf', 'png'):
    fig.savefig(str(OUT) + '/fig_theta_select.' + ext, dpi=600,
                bbox_inches='tight', facecolor='white')
for i, t in enumerate(TH):
    print('theta %-6d n=%d  gen %.4f  IVI %5.0f  [%4.0f, %4.0f]'
          % (t, nn[i], g[i], v[i], vlo[i], vhi[i]))

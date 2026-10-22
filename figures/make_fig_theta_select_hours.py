"""Selecting theta from the generation / reliability trade-off.

Companion to make_fig_theta_select.py, on the reliability axis instead of the
violation axis. One point per theta: annual generation against hours below the
reliability target (below_eps_hrs), DDU, P20, eps = 0.05, Bonferroni.

The good corner is bottom RIGHT: more energy, fewer hours below target.

WHY THIS AXIS BEHAVES AND IVI DOES NOT. Both are read off the same runs, but
hours-below-target counts events and IVI integrates depth, so one harsh inflow
draw moves IVI far more than it moves the count. On the balanced theta = 200
to 1400 set (six paired scenarios) a two-way repeated-measures fit gives:

    hours   residual sd 13.4 h    F(6,30) = 53.4    p = 1.1e-14
    IVI     residual sd 90.6 kaf  F(6,30) =  0.89   p = 0.51

So the theta effect is resolved here and is not resolved on IVI at n = 6.

BAND. Shaded band is the +/- 11.2 h paired 95% confidence interval on a theta
mean, from the residual above (t_0.975,30 * 13.4 / sqrt(6)). It is drawn only
across theta = 200 to 1400, the range where all six scenarios ran. theta = 100,
1600, 3200 and 5000 are at n = 3 and are drawn as open markers; they are not in
the variance estimate and their means are less certain than the band implies.

Output: figures_out/fig_theta_select_hours.{pdf,png}
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
BAL = [200, 400, 600, 800, 1000, 1200, 1400]   # ten scenarios each
PICK = 400
NMAX = 10     # every theta in the refined range now has ten scenarios
CI = 10.0     # h, paired 95% on a theta mean (t_.975,54 * 15.78 / sqrt(10))

g, h, nn = [], [], []
for t in TH:
    d = pd.read_csv(B + 'BPA_TH7_P20_TH%07d/summary_BPA_TH7_P20_TH%07d.csv' % (t, t))
    d = d[(d.policy == 'DDU') & (d.scen <= NMAX)]   # keep the panel BALANCED
    g.append(d.energy_MWh.mean() / 1e6)
    h.append(d.below_eps_hrs.mean())
    nn.append(len(d))
g, h, nn = map(np.array, (g, h, nn))

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

# The path is ordered by theta, NOT by generation: theta = 600 sits slightly
# to the RIGHT of 400 (26.6081 vs 26.6070 TWh), so the trajectory is not a
# function of the x axis and interpolating over x invents a vertical spike
# between them. Straight segments in theta order, and the paired 95% interval
# as bars rather than a band, for the same reason.
ax.plot(g, h, '-', color=C, lw=1.0, alpha=.55, zorder=2)

six, three = nn == NMAX, nn == 3
ax.errorbar(g[six], h[six], yerr=CI, fmt='none', ecolor=C, elinewidth=.9,
            capsize=2.2, capthick=.9, alpha=.85, zorder=3)
ax.plot(g[six], h[six], 'o', ms=4.2, mfc=C, mec=C, mew=1.0, zorder=4,
        label='$n=10$ scenarios, $\\pm$95% CI')
ax.plot(g[three], h[three], 'o', ms=4.2, mfc='white', mec=C, mew=1.1, zorder=4,
        label='$n=3$ scenarios')

k = TH.index(PICK)
ax.plot(g[k], h[k], 'o', ms=6.4, mfc=CH, mec=CH, zorder=5)

OFF = {100: (-2, 9), 200: (6, 4), 400: (-26, -4), 600: (7, -2), 800: (6, -2),
       1000: (-10, -15), 1200: (1, 9), 1400: (-38, -4), 1600: (-4, 9),
       3200: (6, -3), 5000: (-4, 9)}
for i, t in enumerate(TH):
    ax.annotate(r'$\theta=%d$' % t if t in (100, 5000) else r'$%d$' % t,
                (g[i], h[i]), textcoords='offset points', xytext=OFF[t],
                fontsize=7.0, color=CH if t == PICK else INK2)

ax.margins(y=.16, x=.12)
ax.set_xlabel('Annual generation (TWh)')
ax.set_ylabel('Hours below reliability target (h/yr)')
ax.grid(color=RULE, lw=.5)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc='upper left', handletextpad=.4, borderaxespad=.2)

fig.tight_layout(pad=.25)
for ext in ('pdf', 'png'):
    fig.savefig(str(OUT) + '/fig_theta_select_hours.' + ext, dpi=600,
                bbox_inches='tight', facecolor='white')
for i, t in enumerate(TH):
    print('theta %-6d n=%d  gen %.4f TWh  hours %6.1f' % (t, nn[i], g[i], h[i]))

"""Benchmarking: every method lands on one trade-off curve.

Ning's ask: show the chance-constrained solution as a curve over epsilon, with
the robust and scenario solutions marked on the same axes, so a reliability
requirement can be stated and the epsilon that meets it read off.

Axes: annual generation against integrated violation. Up and to the LEFT is
better. Markers only -- the methods are not a sequence, and the reliability
requirement is epsilon itself rather than a line drawn at some violation level,
so neither a connecting line nor a target rule belongs on the panel. All runs are theta = 400, P20 floor, DDU kappa = 1, full water year,
four paired scenarios; the deterministic benchmark has no uncertainty model and
so is the DIU arm by construction.

    chance-constrained   BON, six epsilons (BPA_HL2_TH0000400)
    robust               six budgets, Gamma = 1.0 ... 3.0  (BENCH_TH0000400_E050)
    scenario             six designs, S in {20, 100, 1000}, r in {0, ~5%}
    deterministic        no margin at all

The point of the figure is that these are not four different frontiers. They
are four ways of choosing a point on the same one, because each reduces to a
uniform margin z on the accumulated standard deviation and differs only in how
z is set. What distinguishes the chance constraint is that its knob is the
reliability target itself, so a requirement can be stated and met directly
rather than tuned into.

Output: figures_out/fig_benchmark_eps.{pdf,png}
"""

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures_out'
OUT.mkdir(exist_ok=True)

import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

U = str(ROOT) + '/'
EPS = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25]

# ---- chance-constrained curve
cc = pd.read_csv(U + 'results/runs/BPA_HL2_TH0000400/'
                     'summary_BPA_HL2_TH0000400.csv')
cc = cc[(cc.method == 'BON') & (cc.policy == 'DDU')]
ce, cv = [], []
for e in EPS:
    r = cc[np.isclose(cc.eps, e)]
    ce.append(r.energy_MWh.mean() / 1e3)      # GWh
    cv.append(r.lower_ivi_kaf.mean())
ce, cv = np.array(ce), np.array(cv)

# ---- benchmarks
b = pd.read_csv(U + 'results/runs/benchmarks.csv')


def grp(rule, policy):
    d = b[(b.rule == rule) & (b.policy == policy)]
    out = []
    for (val, nd), s in d.groupby(['val', 'ndisc']):
        out.append((val, nd, s.energy_MWh.mean() / 1e3, s.lower_ivi_kaf.mean()))
    return sorted(out)


rob = grp('robust', 'DDU')
scn = grp('scenario', 'DDU')
det = b[b.rule == 'none']
de, dv = det.energy_MWh.mean() / 1e3, det.lower_ivi_kaf.mean()

CC, ROB, SCN, DET = '#3f7d2a', '#2a78d6', '#cc6600', '#5a5a58'
INK2, RULE = '#555', '#dcdcda'
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 8,
    'axes.labelsize': 8, 'legend.fontsize': 7.0,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'axes.edgecolor': '#666', 'axes.linewidth': .6,
    'xtick.color': '#666', 'ytick.color': '#666', 'axes.labelcolor': 'black',
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white'})

fig, ax = plt.subplots(figsize=(3.5, 3.0))

# Everything that identifies a series lives in the legend. Nothing is written
# onto the panel itself; parameter values go in the caption, in the order the
# markers run.
ax.plot(cv, ce, 'o', ms=4.2, mfc='white', mec=CC, mew=1.2, zorder=4,
        label=r'Chance-constrained, $\epsilon = 0.01$–$0.25$')

ax.plot([r[3] for r in rob], [r[2] for r in rob], 's', ms=3.8, mfc='white',
        mec=ROB, mew=1.1, zorder=4, label=r'Robust, $\Gamma = 1$–$3$')

ax.plot([r[3] for r in scn], [r[2] for r in scn], '^', ms=4.0, mfc='white',
        mec=SCN, mew=1.1, zorder=4, label=r'Scenario, $S = 20$–$1000$')

ax.plot([dv], [de], 'D', ms=4.2, mfc='white', mec=DET, mew=1.1, zorder=4,
        label='Deterministic')

ax.set_xlabel('Integrated violation index (kaf)')
ax.set_ylabel('Annual generation (GWh)')
ax.grid(color=RULE, lw=.5)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc='upper right', handletextpad=.5, fontsize=6.8,
          borderaxespad=.2, labelspacing=.35)

fig.tight_layout(pad=.25)
for ext in ('pdf', 'png'):
    fig.savefig(str(OUT) + '/fig_benchmark_eps.' + ext, dpi=600,
                bbox_inches='tight', facecolor='white')

print('chance-constrained (BON, DDU):')
for e, x, y in zip(EPS, cv, ce):
    print('   eps %.2f   IVI %7.1f   gen %9.0f GWh' % (e, x, y))
print('robust:')
for val, nd, y, x in rob:
    print('   Gamma %-6.3f  IVI %7.1f   gen %9.0f GWh' % (val, x, y))
print('scenario:')
for val, nd, y, x in scn:
    print('   S %-6.0f r %-4.0f IVI %7.1f   gen %9.0f GWh' % (val, nd, x, y))
print('deterministic:      IVI %7.1f   gen %9.0f GWh' % (dv, de))

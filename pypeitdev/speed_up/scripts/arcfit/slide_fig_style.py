"""Shared matplotlib style for the PR-C slide figures (slide_fig_*.py).

Every slide figure is drawn at the exact size (inches) at which it is placed on
the 13.33 x 7.5 in slide, so a 20 pt font in matplotlib is 20 pt on the slide.
No text in these figures is smaller than 20 pt (user rule for the deck).
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEV = os.environ.get('PYPEIT_DEV', '/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite').rstrip('/')
OUTDIR = f'{DEV}/pypeitdev/speed_up/Slides/figures'
INK, MUTED, C_LEG, C_NEW, WIN = '#1b2430', '#6b6a66', '#2a78d6', '#eb6834', '#e9e8e4'
FS = 20


def setup():
    plt.rcParams.update({'font.size': FS, 'axes.titlesize': FS, 'axes.labelsize': FS,
                         'xtick.labelsize': FS, 'ytick.labelsize': FS, 'legend.fontsize': FS,
                         'font.family': 'DejaVu Sans',
                         'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'text.color': INK,
                         'xtick.color': MUTED, 'ytick.color': MUTED,
                         'axes.spines.top': False, 'axes.spines.right': False})


def save(fig, name):
    os.makedirs(OUTDIR, exist_ok=True)
    out = f'{OUTDIR}/{name}'
    # no bbox_inches='tight': keep the figure at its placed size so font sizes are exact
    fig.savefig(out, dpi=200, facecolor='white')
    print('wrote', out)

"""Supplementary figures for the DMG spatial preprint (same style as main set).

figS1: DMG vs GBM patient-level abundance comparison (no significant differences).
figS2: state verification — (a) cross-validated marker AUROC per state per reference
       cohort; (b) marker-score vs deconvolution concordance per section.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': ['Liberation Sans', 'Arimo', 'DejaVu Sans'],
    'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
    'svg.fonttype': 'none', 'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
})
import matplotlib.pyplot as plt

FIG = '/mnt/results/dmg_spatial/figures'
AUD = '/mnt/shared-workspace/audit'
os.makedirs(FIG, exist_ok=True)

STATE_COLORS = {'OPC-like': '#0279EE', 'AC-like': '#FF9400', 'OC-like': '#75A025',
                'MES-like': '#D62728', 'Cycling': '#9467BD', 'Myeloid': '#8C564B',
                'Endothelial': '#17BECF', 'T cell': '#E377C2', 'Oligodendrocyte': '#7F7F7F',
                'Neuron': '#BCBD22', 'Astrocyte': '#FF7F0E', 'Ependymal': '#8C9E5E',
                'Other glia': '#CCCCCC'}


def savefig(fig, name):
    for ext in ('png', 'svg'):
        fig.savefig(f'{FIG}/{name}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('saved', name, flush=True)


def tag(ax, letter):
    ax.text(-0.04, 1.06, letter, transform=ax.transAxes, fontsize=9,
            fontweight='bold', va='top', ha='left')


def figS1():
    d = pd.read_csv(f'{AUD}/dmg_vs_gbm_patient_level.csv')
    d = d.sort_values('DMG_patient_mean').reset_index(drop=True)
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    for i, r in d.iterrows():
        ax.plot([r['DMG_patient_mean'], r['GBM_patient_mean']], [i, i],
                color='#BBBBBB', lw=1.0, zorder=1)
    ax.scatter(d['DMG_patient_mean'], y, s=22, color='#0279EE', zorder=2, label='DMG (n=5 patients)')
    ax.scatter(d['GBM_patient_mean'], y, s=22, color='#FF9400', zorder=2, label='GBM (n=5 patients)')
    ax.set_yticks(y)
    ax.set_yticklabels(d['population'], fontsize=6)
    ax.set_xlabel('Patient-mean abundance (fraction of spots)')
    ax.legend(frameon=False, loc='lower right')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    savefig(fig, 'figS1_dmg_vs_gbm_abundance')


def figS2():
    cv = pd.read_csv(f'{AUD}/cv_auroc_states.csv')
    rob = pd.read_csv(f'{AUD}/marker_vs_c2l_robustness.csv')
    states = ['OPC-like', 'AC-like', 'OC-like', 'MES-like', 'Cycling']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 2.9))
    x = np.arange(len(states))
    coh = {'Liu2022': ('Liu (Smart-seq2)', '#333333'), 'Jessa2022': ('Jessa (10x/snRNA)', '#0279EE')}
    for k, (c, (lbl, col)) in enumerate(coh.items()):
        v = [cv[(cv['cohort'] == c) & (cv['state'] == s)]['cv_AUROC_marker_logreg'].iloc[0] for s in states]
        ax1.bar(x + (k - 0.5) * 0.36, v, width=0.36, color=col, label=lbl)
    ax1.set_xticks(x)
    ax1.set_xticklabels(states, rotation=30, ha='right')
    ax1.set_ylabel('Cross-validated AUROC')
    ax1.set_ylim(0.5, 1.0)
    ax1.axhline(0.5, color='#999999', lw=0.6, ls='--')
    ax1.legend(frameon=False)
    for s in ('top', 'right'):
        ax1.spines[s].set_visible(False)
    tag(ax1, 'a')

    data = [rob[rob['population'] == s]['spearman'].values for s in states]
    bp = ax2.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                     medianprops={'color': 'black', 'lw': 0.8},
                     whiskerprops={'lw': 0.6}, capprops={'lw': 0.6})
    for patch, s in zip(bp['boxes'], states):
        patch.set_facecolor(STATE_COLORS[s]); patch.set_alpha(0.75); patch.set_linewidth(0.6)
    for i, vals in enumerate(data):
        ax2.scatter(np.full(len(vals), i + 1), vals, s=4, color='#333333', zorder=3)
    ax2.set_xticks(range(1, len(states) + 1))
    ax2.set_xticklabels(states, rotation=30, ha='right')
    ax2.set_ylabel('Marker vs deconvolution (Spearman rho)')
    for s in ('top', 'right'):
        ax2.spines[s].set_visible(False)
    tag(ax2, 'b')
    fig.subplots_adjust(wspace=0.35)
    savefig(fig, 'figS2_state_verification')


if __name__ == '__main__':
    figS1()
    figS2()

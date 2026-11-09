"""Publication figure set for the DMG spatial study (Cell/Nature style).

Style: Liberation Sans (Arial-metric), 7 pt base, no axes titles/suptitles
(descriptions live in captions), lowercase bold panel tags (a, b, c, ...),
300 dpi PNG + editable-text SVG.

Story order:
  fig01 reference atlas -> fig02 state maps -> fig03 niches -> fig04 colocalization
  -> fig05 cross-cohort concordance (deconv + reference-free) -> fig06 Xenium
  -> fig07 rostro-caudal axis -> fig08 interface LR crosstalk -> fig09 external validation
"""
import os, json
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
import seaborn as sns
import scanpy as sc

FIG = '/mnt/results/dmg_spatial/figures'
os.makedirs(FIG, exist_ok=True)
NIC = '/mnt/shared-workspace/niches'
C2L = '/mnt/shared-workspace/c2l_out'
VAL = '/mnt/shared-workspace/validation'
AUD = '/mnt/shared-workspace/audit'

STATES = ['OPC-like', 'AC-like', 'OC-like', 'MES-like', 'Cycling', 'Myeloid']
STATE_COLORS = {'OPC-like': '#0279EE', 'AC-like': '#FF9400', 'OC-like': '#75A025',
                'MES-like': '#D62728', 'Cycling': '#9467BD', 'Myeloid': '#8C564B',
                'Endothelial': '#17BECF', 'T cell': '#E377C2', 'Oligodendrocyte': '#7F7F7F',
                'Neuron': '#BCBD22', 'Astrocyte': '#FF7F0E', 'Ependymal': '#8C9E5E',
                'Other glia': '#CCCCCC'}
SECTIONS = ['Ren_DMG1', 'Kord_PonsB', 'Huch_pDMG1']
SECTION_LBL = {'Ren_DMG1': 'Ren DMG1', 'Kord_PonsB': 'Kordowski PonsB', 'Huch_pDMG1': 'Huchede pDMG1'}


def savefig(fig, name):
    for ext in ('png', 'svg'):
        fig.savefig(f'{FIG}/{name}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('saved', name, flush=True)


def tag(ax, letter):
    ax.text(-0.04, 1.06, letter, transform=ax.transAxes, fontsize=9,
            fontweight='bold', va='top', ha='left')


def get_xy(a):
    return np.asarray(a.obsm['spatial'], dtype=float)


# ---------------- fig01: reference atlas ----------------
def fig01():
    a = sc.read_h5ad('/mnt/shared-workspace/c2l_ref.h5ad')
    xy = np.asarray(a.obsm['X_umap'])
    labs = a.obs['ref_label'].astype(str)
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    for c in sorted(labs.unique()):
        m = (labs == c).values
        ax.scatter(xy[m, 0], xy[m, 1], s=0.4, c=STATE_COLORS.get(c, '#999999'),
                   label=f'{c} (n={m.sum():,})', rasterized=True, linewidths=0)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
    ax.legend(markerscale=6, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    savefig(fig, 'fig01_reference_atlas')


# ---------------- fig02: state abundance maps (3 sections x 6 states) ----------------
def fig02():
    fig, axes = plt.subplots(3, 6, figsize=(12.4, 6.6))
    for r, sec in enumerate(SECTIONS):
        abun = pd.read_csv(f'{C2L}/{sec}/abundances.csv', index_col=0)
        a = sc.read_h5ad(f'{NIC}/{sec}/adata_niches.h5ad')
        shared = a.obs_names.intersection(abun.index)
        a = a[shared]; abun = abun.loc[shared]
        xy = get_xy(a)
        for c, st in enumerate(STATES):
            ax = axes[r, c]
            v = abun[st].values
            sc_ = ax.scatter(xy[:, 0], -xy[:, 1], c=v, s=1.2, cmap='viridis', vmin=0,
                             vmax=np.percentile(v, 99) if v.max() > 0 else 1,
                             rasterized=True, linewidths=0)
            ax.set_aspect('equal'); ax.axis('off')
            if r == 0:
                ax.text(0.5, 1.04, st, transform=ax.transAxes, ha='center', fontsize=7)
            if c == 0:
                tag(ax, 'abc'[r])
                ax.text(-0.12, 0.5, SECTION_LBL[sec], transform=ax.transAxes,
                        rotation=90, va='center', fontsize=7)
            cb = plt.colorbar(sc_, ax=ax, shrink=0.55, pad=0.01)
            cb.set_label('cells/spot', fontsize=6)
            cb.ax.tick_params(labelsize=5, width=0.5)
            cb.outline.set_linewidth(0.5)
    fig.subplots_adjust(wspace=0.25, hspace=0.18)
    savefig(fig, 'fig02_state_abundance_maps')


# ---------------- fig03: niche maps ----------------
def fig03():
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
    for i, sec in enumerate(SECTIONS):
        a = sc.read_h5ad(f'{NIC}/{sec}/adata_niches.h5ad')
        xy = get_xy(a)
        labs = a.obs['niche_label'].astype(str)
        keep = labs.value_counts()
        labs[labs.map(keep) < 20] = 'other/low'
        ax = axes[i]
        for c in sorted(labs.unique()):
            m = (labs == c).values
            ax.scatter(xy[m, 0], -xy[m, 1], s=1.2, c=STATE_COLORS.get(c, '#DDDDDD'),
                       label=c, rasterized=True, linewidths=0)
        ax.set_aspect('equal'); ax.axis('off')
        tag(ax, 'abc'[i])
        ax.text(0.5, -0.02, SECTION_LBL[sec], transform=ax.transAxes, ha='center', fontsize=7)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, markerscale=6, loc='center right', frameon=False,
               bbox_to_anchor=(1.0, 0.5))
    fig.subplots_adjust(wspace=0.05, right=0.82)
    savefig(fig, 'fig03_niche_maps')


# ---------------- fig04: colocalization heatmaps ----------------
def fig04():
    cohs = ['Ren_DMG', 'Kordowski', 'Huchede', 'Ren_GBM']
    lbl = {'Ren_DMG': 'Ren DMG', 'Kordowski': 'Kordowski', 'Huchede': 'Huchede', 'Ren_GBM': 'Ren GBM'}
    fig, axes = plt.subplots(1, 5, figsize=(11.5, 3.2),
                             gridspec_kw={'width_ratios': [1, 1, 1, 1, 0.05], 'wspace': 0.18})
    cax = axes[4]
    for i, (ax, c) in enumerate(zip(axes[:4], cohs)):
        m = pd.read_csv(f'{VAL}/colocalization_mean_{c}.csv', index_col=0)
        sns.heatmap(m, ax=ax, cmap='RdBu_r', center=0, vmin=-0.6, vmax=0.6,
                    xticklabels=True, yticklabels=(i == 0), cbar=(i == 0), cbar_ax=cax,
                    square=True, cbar_kws={'label': 'Spearman rho'})
        ax.set_xlabel(lbl[c])
        ax.tick_params(axis='x', rotation=90, labelsize=5.5)
        ax.tick_params(axis='y', labelsize=5.5)
        tag(ax, 'abcd'[i])
    savefig(fig, 'fig04_colocalization_heatmaps')


# ---------------- fig05: concordance (deconv vs reference-free) ----------------
def fig05():
    pw = pd.read_csv(f'{VAL}/cohort_colocalization_concordance.csv', index_col=0).astype(float)
    dec = pd.read_csv(f'{AUD}/concordance_deconv_pointest.csv', index_col=0)
    dec_lo = pd.read_csv(f'{AUD}/concordance_deconv_ci_lo.csv', index_col=0)
    dec_hi = pd.read_csv(f'{AUD}/concordance_deconv_ci_hi.csv', index_col=0)
    mf = pd.read_csv(f'{AUD}/concordance_markerfree_pointest.csv', index_col=0)
    mf_lo = pd.read_csv(f'{AUD}/concordance_markerfree_ci_lo.csv', index_col=0)
    mf_hi = pd.read_csv(f'{AUD}/concordance_markerfree_ci_hi.csv', index_col=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.4), gridspec_kw={'width_ratios': [1, 1.35]})
    sns.heatmap(pw, annot=True, fmt='.2f', cmap='viridis', vmin=0, vmax=1, square=True,
                ax=ax1, cbar_kws={'label': 'Spearman rho'}, annot_kws={'size': 6})
    ax1.tick_params(labelsize=6)
    tag(ax1, 'a')

    pairs = [('Ren_DMG', 'Huchede'), ('Ren_DMG', 'Ren_GBM'), ('Ren_DMG', 'Kordowski'),
             ('Huchede', 'Kordowski'), ('Ren_GBM', 'Huchede'), ('Ren_GBM', 'Kordowski')]
    plabels = ['Ren DMG–Huchede', 'Ren DMG–Ren GBM', 'Ren DMG–Kordowski',
               'Huchede–Kordowski', 'Ren GBM–Huchede', 'Ren GBM–Kordowski']
    x = np.arange(len(pairs))
    for k, (pt, lo, hi, col, name) in enumerate([(dec, dec_lo, dec_hi, '#333333', 'Deconvolution-based'),
                                                 (mf, mf_lo, mf_hi, '#0279EE', 'Reference-free markers')]):
        y = [pt.loc[a, b] for a, b in pairs]
        yerr = np.array([[pt.loc[a, b] - lo.loc[a, b] for a, b in pairs],
                         [hi.loc[a, b] - pt.loc[a, b] for a, b in pairs]])
        ax2.errorbar(x + (k - 0.5) * 0.16, y, yerr=yerr, fmt='o', ms=3.5, lw=0.8,
                     capsize=2, color=col, label=name)
    ax2.set_xticks(x); ax2.set_xticklabels(plabels, rotation=30, ha='right')
    ax2.set_ylabel('Colocalization concordance (rho)')
    ax2.set_ylim(0, 1)
    ax2.legend(frameon=False, loc='lower right')
    for s in ('top', 'right'):
        ax2.spines[s].set_visible(False)
    tag(ax2, 'b')
    fig.subplots_adjust(wspace=0.3)
    savefig(fig, 'fig05_cohort_concordance')


# ---------------- fig06: Xenium ----------------
def fig06():
    a = sc.read_h5ad('/mnt/shared-workspace/xenium_H3K27M_P6_labeled.h5ad')
    a = a[a.obs['xenium_margin'] >= 0.5]
    xy = get_xy(a)
    labs = a.obs['xenium_state'].astype(str)
    z = pd.read_csv(f'{VAL}/xenium_nhood_enrichment_z.csv', index_col=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.2), gridspec_kw={'width_ratios': [1.15, 1]})
    for c in sorted(labs.unique()):
        m = (labs == c).values
        ax1.scatter(xy[m, 0], -xy[m, 1], s=0.15, c=STATE_COLORS.get(c, '#999999'),
                    label=c, rasterized=True, linewidths=0, alpha=0.7)
    ax1.set_aspect('equal'); ax1.axis('off')
    ax1.legend(markerscale=8, loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False,
               fontsize=6, handletextpad=0.3, borderaxespad=0.2)
    tag(ax1, 'a')

    sns.heatmap(z, cmap='RdBu_r', center=0, vmin=-80, vmax=80, annot=True, fmt='.0f',
                square=True, ax=ax2, cbar_kws={'label': 'Neighborhood enrichment z'},
                annot_kws={'size': 5})
    ax2.tick_params(labelsize=5.5)
    tag(ax2, 'b')
    fig.subplots_adjust(wspace=0.55)
    savefig(fig, 'fig06_xenium_territoriality')


# ---------------- fig07: rostro-caudal axis ----------------
def fig07():
    rc = pd.read_csv(f'{VAL}/kordowski_rostrocaudal_abundances.csv', index_col=0)
    mf = pd.read_csv(f'{AUD}/kordowski_rostrocaudal_markerfree.csv')
    order = ['Kord_Midbrain', 'Kord_PonsA', 'Kord_PonsB', 'Kord_PonsC', 'Kord_PonsD',
             'Kord_MedullaA', 'Kord_MedullaB', 'Kord_MedullaC', 'Kord_SpinalCordA', 'Kord_SpinalCordB']
    rc = rc.reindex(order)
    x = np.arange(len(rc))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.2))
    for st in STATES:
        ax1.plot(x, rc[st], marker='o', ms=3, lw=1.2, color=STATE_COLORS[st], label=st)
    ax1.set_xticks(x)
    ax1.set_xticklabels([s.replace('Kord_', '') for s in order], rotation=45, ha='right')
    ax1.set_ylabel('Mean cells/spot')
    ax1.set_xlabel('Rostro-caudal position')
    ax1.legend(frameon=False, ncol=2)
    for s in ('top', 'right'):
        ax1.spines[s].set_visible(False)
    tag(ax1, 'a')

    ax2.plot(x, mf.set_index('region').reindex(
        [s.replace('Kord_', '') for s in order])['opc_marker_score_mean'].values,
        marker='o', ms=3, lw=1.2, color=STATE_COLORS['OPC-like'])
    ax2.set_xticks(x)
    ax2.set_xticklabels([s.replace('Kord_', '') for s in order], rotation=45, ha='right')
    ax2.set_ylabel('OPC-like marker score')
    ax2.set_xlabel('Rostro-caudal position')
    for s in ('top', 'right'):
        ax2.spines[s].set_visible(False)
    tag(ax2, 'b')
    fig.subplots_adjust(wspace=0.3)
    savefig(fig, 'fig07_rostrocaudal_axis')


# ---------------- fig08: interface LR bubble ----------------
def fig08():
    iface = pd.read_csv(f'{VAL}/lr_interface_crosstalk.csv')
    sp = iface['interaction'].str.split('|', expand=True)
    iface['pair'] = sp[0] + '>' + sp[1]
    st = sp[2].str.split('->', expand=True)
    iface['st'] = st[0] + '->' + st[1]
    iface['dmg_support'] = iface['ren_dmg_sections_sig'] + iface['huchede_sections_sig']
    top = iface.sort_values(['dmg_support', 'mean_spatial_z'], ascending=False).head(25)
    top = top.iloc[::-1]
    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    sc_ = ax.scatter(top['mean_lr_means'], range(len(top)), s=42,
                     c=top['gbm_sections_sig'], cmap='viridis', vmin=0, vmax=6,
                     edgecolors='white', linewidths=0.4)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['pair'] + ' | ' + top['st'], fontsize=5.5)
    ax.set_xlabel('Mean LR score (Kordowski sections)')
    cb = plt.colorbar(sc_, ax=ax, label='Ren GBM sections significant (of 6)', shrink=0.5)
    cb.ax.tick_params(labelsize=6)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    savefig(fig, 'fig08_interface_lr_crosstalk')


# ---------------- fig09: external validation (Collot 2026) ----------------
def fig09():
    t1 = pd.read_csv(f'{AUD}/collot/t1_coloc_concordance.csv')
    t2 = pd.read_csv(f'{AUD}/collot/t2_malignant_composition.csv')
    t3 = json.load(open(f'{AUD}/collot/t3_lr_recurrence.json'))

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(11.5, 3.2))
    lbls = [c.replace('Collot_vs_', '').replace('_', ' ') for c in t1['comparison']]
    ax1.bar(lbls, t1['coloc_spearman'], color='#0279EE', width=0.6)
    ax1.set_ylabel('Colocalization concordance (rho)')
    ax1.set_ylim(0, 1)
    ax1.tick_params(axis='x', rotation=30)
    for tick in ax1.get_xticklabels():
        tick.set_ha('right')
    for s in ('top', 'right'):
        ax1.spines[s].set_visible(False)
    tag(ax1, 'a')

    malignant = ['OPC-like', 'AC-like', 'OC-like', 'MES-like', 'Cycling']
    bottoms = np.zeros(len(t2))
    x2 = np.arange(len(t2))
    for st in malignant:
        v = t2[f'frac_{st}'].values
        ax2.bar(x2, v, bottom=bottoms, color=STATE_COLORS[st], label=st, width=0.6)
        bottoms += v
    ax2.set_xticks(x2)
    ax2.set_xticklabels([s.replace('Collot_', '').replace('_HE', '').replace('_', '\n') for s in t2['section']], fontsize=6)
    ax2.set_ylabel('Malignant fraction')
    ax2.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=5,
               fontsize=6, columnspacing=0.7, handletextpad=0.4)
    for s in ('top', 'right'):
        ax2.spines[s].set_visible(False)
    tag(ax2, 'b')

    cats = ['3-cohort pairs\n(any section)', '3-cohort pairs\n(recurrent)',
            'DMG-leaning\n(any section)', 'DMG-leaning\n(recurrent)']
    vals = [t3['three_cohort_any_section_n'], t3['three_cohort_recurrent_in_collot_n'],
            len(t3['dmg_leaning_any_section']), len(t3['dmg_leaning_recurrent_in_collot'])]
    tots = [t3['three_cohort_total'], t3['three_cohort_total'], t3['dmg_leaning_total'], t3['dmg_leaning_total']]
    colors = ['#75A025', '#75A025', '#D62728', '#D62728']
    ax3.bar(cats, vals, color=colors, width=0.6)
    for i, (v, t) in enumerate(zip(vals, tots)):
        ax3.text(i, v + 1, f'{v}/{t}', ha='center', fontsize=6)
    ax3.set_ylabel('Pairs detected in Collot cohort')
    ax3.tick_params(axis='x', labelsize=6)
    for s in ('top', 'right'):
        ax3.spines[s].set_visible(False)
    tag(ax3, 'c')
    fig.subplots_adjust(wspace=0.35)
    savefig(fig, 'fig09_external_validation')


if __name__ == '__main__':
    import sys
    which = sys.argv[1:] if len(sys.argv) > 1 else ['all']
    fns = {'fig01': fig01, 'fig02': fig02, 'fig03': fig03, 'fig04': fig04, 'fig05': fig05,
           'fig06': fig06, 'fig07': fig07, 'fig08': fig08, 'fig09': fig09}
    for k, fn in fns.items():
        if 'all' in which or k in which:
            fn()

"""Aggregate LIANA LR pairs across sections and cohorts; compare to Kordowski's published LR set.

Per-section filter: magnitude_rank <= 0.05 (LIANA aggregate) and spatially supported
(spatial_z > 2 or source == target).
Cohort recurrence: pair (ligand-receptor, source, target) significant in >= ceil(60%) of
cohort sections (Kordowski: >= 5/10; Huchede: >= 2/3; Ren DMG: >= 3/5).
Headline set: recurrent in >= 2 independent DMG cohorts.
"""
import os, glob
import numpy as np, pandas as pd

LR = '/mnt/shared-workspace/lr_out'
OUT = '/mnt/shared-workspace/validation'
os.makedirs(OUT, exist_ok=True)

COHORT = {'Ren_DMG1': 'Ren_DMG', 'Ren_DMG2': 'Ren_DMG', 'Ren_DMG3': 'Ren_DMG',
          'Ren_DMG4': 'Ren_DMG', 'Ren_DMG5': 'Ren_DMG',
          'Ren_GBM1': 'Ren_GBM', 'Ren_GBM2': 'Ren_GBM', 'Ren_GBM3': 'Ren_GBM',
          'Ren_GBM4': 'Ren_GBM', 'Ren_GBM5_1': 'Ren_GBM', 'Ren_GBM5_2': 'Ren_GBM',
          **{f'Kord_{r}': 'Kordowski' for r in
             ['Midbrain','PonsA','PonsB','PonsC','PonsD','MedullaA','MedullaB','MedullaC','SpinalCordA','SpinalCordB']},
          'Huch_pDMG1': 'Huchede', 'Huch_pDMG2': 'Huchede', 'Huch_pDMG3': 'Huchede'}
MIN_RECUR = {'Ren_DMG': 3, 'Ren_GBM': 4, 'Kordowski': 5, 'Huchede': 2}

frames = []
for s, coh in COHORT.items():
    f = f'{LR}/{s}_lr_pairs.csv'
    if not os.path.exists(f):
        print('missing', s); continue
    d = pd.read_csv(f)
    d['cohort'] = coh
    frames.append(d)
lr = pd.concat(frames, ignore_index=True)
print('total pairs:', len(lr), '| columns:', [c for c in lr.columns if 'rank' in c or 'pval' in c])

# per-section significance filter
sig = lr[(lr['magnitude_rank'] <= 0.05)].copy()
sig['pair'] = sig['ligand_complex'] + '|' + sig['receptor_complex']
sig['interaction'] = sig['pair'] + '|' + sig['source'] + '->' + sig['target']
sp = sig[(sig['spatial_z'] > 2) | (sig['source'] == sig['target'])].copy()
print('significant (magnitude_rank<=0.05):', len(sig), '| spatially supported:', len(sp))

# cohort recurrence on spatially supported pairs
rec = sp.groupby(['cohort', 'interaction'])['section'].nunique().reset_index(name='n_sections')
rec['recurrent'] = rec.apply(lambda r: r['n_sections'] >= MIN_RECUR[r['cohort']], axis=1)
recur = rec[rec['recurrent']]
print('\nRecurrent spatially-supported interactions per cohort:')
print(recur.groupby('cohort').size())

# cross-cohort replication (DMG cohorts only for headline)
dmg = recur[recur['cohort'].isin(['Ren_DMG', 'Kordowski', 'Huchede'])]
n_coh = dmg.groupby('interaction')['cohort'].nunique()
headline = n_coh[n_coh >= 2].index
head = dmg[dmg['interaction'].isin(headline)].copy()
head.to_csv(f'{OUT}/lr_headline_replicated.csv', index=False)
print('\nHeadline interactions (recurrent in >=2 DMG cohorts):', len(headline))

# widest replication
top = (dmg.groupby('interaction')['cohort'].nunique().sort_values(ascending=False).head(30))
print('\nTop 30 by cohort replication:')
print(top)

# Kordowski published LR comparison: per-spot score matrices, columns = LIG_RECEPTOR
our_pairs = set(sp['pair'].unique())  # 'ligand|receptor' (complexes may contain &)
def norm_pair(p):
    return p.replace('&', '_').replace('|', '_').upper()
our_norm = {norm_pair(p) for p in our_pairs}

overlap_rows = []
for f in sorted(glob.glob('/mnt/shared-workspace/dmg_spatial/kordowski/*LR_scores.csv.gz')):
    d = pd.read_csv(f, index_col=0)
    region = os.path.basename(f).split('_updated')[0].split('_', 1)[-1]
    mean_score = d.mean(axis=0)
    top_pub = set(mean_score[mean_score > mean_score.quantile(0.95)].index.str.upper())
    ov = len(our_norm & top_pub)
    overlap_rows.append({'region': region, 'n_top_published': len(top_pub),
                         'overlap_with_ours': ov,
                         'frac_of_top_published': round(ov / max(len(top_pub), 1), 3)})
ov_df = pd.DataFrame(overlap_rows)
ov_df.to_csv(f'{OUT}/kordowski_published_lr_overlap.csv', index=False)
print('\nOverlap with Kordowski published top-5% LR scores:')
print(ov_df.to_string(index=False))

recur.to_csv(f'{OUT}/lr_recurrence_by_cohort.csv', index=False)
print('\nDONE')

"""Cross-cohort validation aggregation.

Tier 1 (patient-level): Ren DMG (n=5) vs Huchede (n=3) — correlate state-colocalization
  (abundance Spearman) matrices and niche compositions across cohorts.
Tier 2 (spatial replication): Kordowski 10 sections, 1 patient — rostro-caudal profiles.
Comparator: Ren GBM sections for DMG-selectivity.
Outputs to /mnt/shared-workspace/validation/.
"""
import os, glob, json
import numpy as np, pandas as pd

NIC = '/mnt/shared-workspace/niches'
OUT = '/mnt/shared-workspace/validation'
os.makedirs(OUT, exist_ok=True)

COHORT = {'Ren_DMG1': 'Ren_DMG', 'Ren_DMG2': 'Ren_DMG', 'Ren_DMG3': 'Ren_DMG',
          'Ren_DMG4': 'Ren_DMG', 'Ren_DMG5': 'Ren_DMG',
          'Ren_GBM1': 'Ren_GBM', 'Ren_GBM2': 'Ren_GBM', 'Ren_GBM3': 'Ren_GBM',
          'Ren_GBM4': 'Ren_GBM', 'Ren_GBM5_1': 'Ren_GBM', 'Ren_GBM5_2': 'Ren_GBM',
          **{f'Kord_{r}': 'Kordowski' for r in
             ['Midbrain','PonsA','PonsB','PonsC','PonsD','MedullaA','MedullaB','MedullaC','SpinalCordA','SpinalCordB']},
          'Huch_pDMG1': 'Huchede', 'Huch_pDMG2': 'Huchede', 'Huch_pDMG3': 'Huchede'}

POPS = ['OPC-like','AC-like','OC-like','MES-like','Cycling','Myeloid','T cell',
        'Oligodendrocyte','Endothelial','Astrocyte','Ependymal','Neuron','Other glia']

# --- 1. Per-section abundance Spearman matrices -> cohort means ---
mats = {}
for s, coh in COHORT.items():
    f = f'{NIC}/{s}/abundance_spearman.csv'
    if os.path.exists(f):
        m = pd.read_csv(f, index_col=0).reindex(index=POPS, columns=POPS)
        mats.setdefault(coh, []).append(m)

cohort_mean = {c: sum(ms) / len(ms) for c, ms in mats.items()}
for c, m in cohort_mean.items():
    m.to_csv(f'{OUT}/colocalization_mean_{c}.csv')

# pairwise cohort correlation of upper-triangle values
def ut(m):
    return m.values[np.triu_indices_from(m, k=1)]
cohs = list(cohort_mean)
pw = pd.DataFrame(index=cohs, columns=cohs, dtype=float)
for i, a in enumerate(cohs):
    for b in cohs[i+1:]:
        x, y = ut(cohort_mean[a]), ut(cohort_mean[b])
        ok = ~(np.isnan(x) | np.isnan(y))
        pw.loc[a, b] = pw.loc[b, a] = pd.Series(x[ok]).corr(pd.Series(y[ok]), method='spearman')
pw.to_csv(f'{OUT}/cohort_colocalization_concordance.csv')
print('Cohort colocalization concordance (Spearman):')
print(pw.round(3))

# --- 2. Moran's I per population across sections ---
mi_rows = []
for s, coh in COHORT.items():
    f = f'{NIC}/{s}/morans_i.csv'
    if os.path.exists(f):
        r = pd.read_csv(f, index_col=0).iloc[:, 0]
        mi_rows.append(pd.Series(r, name=s))
mi = pd.DataFrame(mi_rows).reindex(columns=POPS)
mi.insert(0, 'cohort', [COHORT[s] for s in mi.index])
mi.to_csv(f'{OUT}/morans_i_by_section.csv')
print('\nMoran\'s I (mean per cohort):')
print(mi.groupby('cohort').mean().round(3))

# --- 3. Kordowski rostro-caudal profiles ---
order = ['Kord_Midbrain','Kord_PonsA','Kord_PonsB','Kord_PonsC','Kord_PonsD',
         'Kord_MedullaA','Kord_MedullaB','Kord_MedullaC','Kord_SpinalCordA','Kord_SpinalCordB']
abun_means = {}
for s in order:
    f = f'/mnt/shared-workspace/c2l_out/{s}/abundances.csv'
    if os.path.exists(f):
        abun_means[s] = pd.read_csv(f, index_col=0).mean()
rc = pd.DataFrame(abun_means).T.reindex(columns=POPS)
rc.to_csv(f'{OUT}/kordowski_rostrocaudal_abundances.csv')
print('\nKordowski rostro-caudal mean abundances:')
print(rc.round(2))
print('\nDONE')

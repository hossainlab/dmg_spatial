"""Fix spatial support for LR pairs: direct neighborhood enrichment between
dominant-population labels (instead of lossy niche-level mapping).
Reuses cached LIANA results; recomputes spatial_z and re-filters.
"""
import sys, os
import numpy as np, pandas as pd
import scanpy as sc
import squidpy as sq

C2L = '/mnt/shared-workspace/c2l_out'
QC = '/mnt/shared-workspace/qc_domains'
LR = '/mnt/shared-workspace/lr_out'


def fix_section(name):
    lr = pd.read_csv(f'{LR}/{name}_lr_pairs.csv')
    abun = pd.read_csv(f'{C2L}/{name}/abundances.csv', index_col=0)
    a = sc.read_h5ad(f'{QC}/{name}/adata_domains.h5ad')
    shared = a.obs_names.intersection(abun.index)
    a = a[shared].copy()
    a.obs['dominant'] = abun.loc[shared].idxmax(axis=1).values
    vc = a.obs['dominant'].value_counts()
    keep = vc[vc >= 20].index
    a = a[a.obs['dominant'].isin(keep)].copy()
    a.obs['dominant'] = a.obs['dominant'].astype('category')
    if 'spatial_connectivities' not in a.obsp:
        sq.gr.spatial_neighbors(a, coord_type='grid', n_neighs=6)
    sq.gr.nhood_enrichment(a, cluster_key='dominant')
    cats = [str(c) for c in a.obs['dominant'].cat.categories]
    z = pd.DataFrame(a.uns['dominant_nhood_enrichment']['zscore'], index=cats, columns=cats)
    lr['spatial_z'] = [z.loc[s, t] if (s in z.index and t in z.columns) else np.nan
                       for s, t in zip(lr['source'], lr['target'])]
    lr.to_csv(f'{LR}/{name}_lr_pairs.csv', index=False)
    sp = lr[(lr['spatial_z'] > 2) | (lr['source'] == lr['target'])]
    sp.to_csv(f'{LR}/{name}_lr_pairs_spatial.csv', index=False)
    return {'section': name, 'n_pairs': len(lr), 'n_spatial': len(sp),
            'n_spatial_cross': int((sp['source'] != sp['target']).sum())}


if __name__ == '__main__':
    rows = [fix_section(n) for n in sys.argv[1:]]
    print(pd.DataFrame(rows).to_string(index=False))

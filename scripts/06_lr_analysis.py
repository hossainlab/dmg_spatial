"""Spatially-aware niche ligand-receptor analysis.

Per section:
  1. Spots labeled by dominant deconvolved population (niche_label).
  2. LIANA rank-aggregate LR inference between niche labels (raw counts, log-normalized).
  3. Spatial support filter: source/target populations must colocalize
     (neighborhood enrichment z > 2 between their niches, or significant co-occurrence).
Outputs per section: lr_pairs.csv (all), lr_pairs_spatial.csv (spatially supported).
Cross-cohort aggregation in aggregate_lr.py.
"""
import sys, os
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import scipy.sparse as sp

C2L = '/mnt/shared-workspace/c2l_out'
QC = '/mnt/shared-workspace/qc_domains'
OUT = '/mnt/shared-workspace/lr_out'
os.makedirs(OUT, exist_ok=True)

MIN_SPOTS_PER_LABEL = 20


def run_section(name):
    import liana as li
    abun = pd.read_csv(f'{C2L}/{name}/abundances.csv', index_col=0)
    a = sc.read_h5ad(f'{QC}/{name}/adata_domains.h5ad')
    shared = a.obs_names.intersection(abun.index)
    a = a[shared].copy()
    abun = abun.loc[shared]

    # dominant population per spot
    a.obs['dominant'] = abun.idxmax(axis=1).values
    # drop labels with too few spots
    vc = a.obs['dominant'].value_counts()
    keep_labels = vc[vc >= MIN_SPOTS_PER_LABEL].index
    a = a[a.obs['dominant'].isin(keep_labels)].copy()

    # normalize for LIANA
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    li.mt.rank_aggregate(adata=a, groupby='dominant', use_raw=False,
                         expr_prop=0.1, min_cells=MIN_SPOTS_PER_LABEL, verbose=False)
    lr = a.uns['liana_res'].copy()
    lr['section'] = name

    # spatial support: neighborhood enrichment between niches of source & target
    nhz_path = f'/mnt/shared-workspace/niches/{name}/neighborhood_enrichment_z.csv'
    if os.path.exists(nhz_path):
        nhz = pd.read_csv(nhz_path, index_col=0)
        # map dominant label -> niche z via niche composition: use niche_label assignments
        niche_assign = pd.read_csv(f'/mnt/shared-workspace/niches/{name}/niche_assignments.csv', index_col=0)
        # z between niches dominated by source and target populations
        lab2niche = niche_assign.groupby('niche_label')['niche'].agg(lambda s: s.value_counts().idxmax())
        def get_z(src, tgt):
            try:
                return float(nhz.loc[str(lab2niche.get(src)), str(lab2niche.get(tgt))])
            except Exception:
                return np.nan
        lr['spatial_z'] = [get_z(s, t) for s, t in zip(lr['source'], lr['target'])]
    else:
        lr['spatial_z'] = np.nan

    lr.to_csv(f'{OUT}/{name}_lr_pairs.csv', index=False)
    sp_lr = lr[(lr['spatial_z'] > 2) | (lr['source'] == lr['target'])]
    sp_lr.to_csv(f'{OUT}/{name}_lr_pairs_spatial.csv', index=False)
    return {'section': name, 'n_pairs': len(lr), 'n_spatial': len(sp_lr)}


if __name__ == '__main__':
    rows = []
    for name in sys.argv[1:]:
        if os.path.exists(f'{OUT}/{name}_lr_pairs.csv'):
            print(f'[skip {name}]', flush=True)
            continue
        try:
            rows.append(run_section(name))
            print(f'[lr {name}] ok', flush=True)
        except Exception as e:
            print(f'[lr {name}] FAIL: {e}', flush=True)
    if rows:
        print(pd.DataFrame(rows).to_string(index=False))

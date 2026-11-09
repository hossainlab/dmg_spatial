"""Niche definition + colocalization statistics from cell2location abundances.

Per section:
  1. Join c2l q05 abundances with QC'd spatial AnnData (by barcode).
  2. Composition-Leiden niche clustering on neighbor-smoothed abundances.
  3. Colocalization: neighborhood enrichment (squidpy) on niches, Moran's I per
     population, pairwise spatial Spearman of abundance maps.
Outputs per section: niche_assignments.csv, neighborhood_enrichment_z.csv,
morans_i.csv, abundance_spearman.csv, niche_composition.csv, adata_niches.h5ad.
"""
import sys, os, json, shutil
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import scipy.sparse as sp
import squidpy as sq

C2L_OUT = '/mnt/shared-workspace/c2l_out'
QC_OUT = '/mnt/shared-workspace/qc_domains'
NICHE_OUT = '/mnt/shared-workspace/niches'
os.makedirs(NICHE_OUT, exist_ok=True)


def moran_i(x, C):
    x = x - x.mean()
    n = len(x)
    S0 = C.sum()
    den = (x ** 2).sum()
    return (n / S0) * ((C @ x * x).sum() / den) if den > 0 and S0 > 0 else np.nan


def analyze_section(name, smooth_rounds=2):
    abun = pd.read_csv(f'{C2L_OUT}/{name}/abundances.csv', index_col=0)
    a = sc.read_h5ad(f'{QC_OUT}/{name}/adata_domains.h5ad')
    shared = a.obs_names.intersection(abun.index)
    a = a[shared].copy()
    abun = abun.loc[shared]
    for c in abun.columns:
        a.obs[f'abun_{c}'] = abun[c].values

    if 'spatial_connectivities' not in a.obsp:
        sq.gr.spatial_neighbors(a, coord_type='grid', n_neighs=6)
    conn = a.obsp['spatial_connectivities']

    # neighbor-smoothed composition -> niche clustering
    Xs = abun.values.copy()
    for _ in range(smooth_rounds):
        Xs = 0.5 * Xs + 0.5 * (conn @ Xs) / np.maximum(conn.sum(1).A1[:, None], 1)
    an = ad.AnnData(X=np.clip(Xs, 0, None))
    sc.pp.normalize_total(an, target_sum=1)
    sc.pp.log1p(an)
    sc.pp.pca(an, n_comps=min(10, Xs.shape[1] - 1))
    sc.pp.neighbors(an, n_neighbors=15)
    sc.tl.leiden(an, resolution=0.8)
    a.obs['niche'] = an.obs['leiden'].values

    comp = a.obs.groupby('niche')[[f'abun_{c}' for c in abun.columns]].mean()
    comp.columns = [c.replace('abun_', '') for c in comp.columns]
    niche_label = comp.idxmax(axis=1)
    a.obs['niche_label'] = a.obs['niche'].map(niche_label)

    sq.gr.nhood_enrichment(a, cluster_key='niche')
    cats = [str(c) for c in a.obs['niche'].cat.categories]
    nhood_z = pd.DataFrame(a.uns['niche_nhood_enrichment']['zscore'], index=cats, columns=cats)

    morans = {c: moran_i(a.obs[f'abun_{c}'].values, conn) for c in abun.columns
              if a.obs[f'abun_{c}'].std() > 0}
    corr = abun.corr(method='spearman')

    sdir = f'{NICHE_OUT}/{name}'
    os.makedirs(sdir, exist_ok=True)
    os.makedirs('/workspace/niches_tmp', exist_ok=True)
    a.obs[['niche', 'niche_label']].to_csv(f'{sdir}/niche_assignments.csv')
    nhood_z.to_csv(f'{sdir}/neighborhood_enrichment_z.csv')
    pd.Series(morans).to_csv(f'{sdir}/morans_i.csv')
    corr.to_csv(f'{sdir}/abundance_spearman.csv')
    comp.to_csv(f'{sdir}/niche_composition.csv')
    tmp = f'/workspace/niches_tmp/{name}_adata_niches.h5ad'
    a.write(tmp, compression='gzip')
    shutil.copy(tmp, f'{sdir}/adata_niches.h5ad')
    return {'section': name, 'n_niches': a.obs['niche'].nunique(),
            'morans_OPC': round(morans.get('OPC-like', np.nan), 3)}


if __name__ == '__main__':
    rows = []
    for name in sys.argv[1:]:
        if os.path.exists(f'{NICHE_OUT}/{name}/niche_assignments.csv'):
            print(f'[niche {name}] skip (done)', flush=True)
            continue
        try:
            rows.append(analyze_section(name))
            print(f'[niche {name}] ok', flush=True)
        except Exception as e:
            print(f'[niche {name}] FAIL: {e}', flush=True)
            rows.append({'section': name, 'n_niches': None, 'morans_OPC': np.nan})
    if rows:
        print(pd.DataFrame(rows).to_string(index=False))

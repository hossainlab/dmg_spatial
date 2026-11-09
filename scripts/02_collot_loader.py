"""Loader for Collot & Ruiz-Moreno et al. 2026 (Cancer Cell) Visium Space Ranger outputs.

Zenodo 10.5281/zenodo.17922773, ST.tar.gz — newly generated primary human DMG
samples only (per README; Ren/Kordowski/Van den Broek re-analyses NOT included).
"""
import os, glob
import pandas as pd, numpy as np
import scanpy as sc
import scipy.sparse as sp

COLLOT_DIR = '/mnt/shared-workspace/collot/ST'


def list_samples():
    """Return {sample_name: spaceranger_dir} for every human Visium sample."""
    out = {}
    for d in sorted(glob.glob(f'{COLLOT_DIR}/*')):
        if not os.path.isdir(d):
            continue
        # Space Ranger dir must contain filtered_feature_bc_matrix
        if glob.glob(f'{d}/**/filtered_feature_bc_matrix', recursive=True):
            out[os.path.basename(d)] = d
    return out


def load_collot(sample):
    """Load one Collot Space Ranger output -> AnnData with raw counts + spatial."""
    d = list_samples()[sample]
    fbc = glob.glob(f'{d}/**/filtered_feature_bc_matrix', recursive=True)[0]
    a = sc.read_10x_mtx(fbc)
    a.var_names_make_unique()

    # spatial coordinates (tissue_positions_list.csv for SR<3, tissue_positions.csv for SR>=3)
    spdir = glob.glob(f'{d}/**/spatial', recursive=True)[0]
    pos_path = None
    for cand in ('tissue_positions_list.csv', 'tissue_positions.csv'):
        p = os.path.join(spdir, cand)
        if os.path.exists(p):
            pos_path = p
            break
    pos = pd.read_csv(pos_path, header=None if 'list' in pos_path else 0)
    if 'barcode' not in pos.columns:
        pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col',
                       'pxl_row_in_fullres', 'pxl_col_in_fullres']
    pos = pos.set_index('barcode')
    a.obs = a.obs.join(pos, how='left')
    a.obs['in_tissue'] = a.obs['in_tissue'].fillna(0).astype(int)
    a = a[a.obs['in_tissue'] == 1].copy()
    a.obsm['spatial'] = a.obs[['pxl_col_in_fullres', 'pxl_row_in_fullres']].values
    a.obs['array_row'] = a.obs['array_row'].astype(int)
    a.obs['array_col'] = a.obs['array_col'].astype(int)
    a.obs['sample'] = sample
    a.obs['cohort_spatial'] = 'Collot2026'
    a.obs['group'] = 'DMG'
    a.X = sp.csr_matrix(a.X)
    return a


if __name__ == '__main__':
    s = list_samples()
    print(f'{len(s)} samples:')
    for k, v in s.items():
        print(f'  {k}: {v}')

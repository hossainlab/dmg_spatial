"""Loaders for the DMG spatial datasets (Ren / Kordowski / Huchede Visium, Xenium)."""
import os, glob, json
import pandas as pd, numpy as np
import scanpy as sc, anndata as ad
import scipy.sparse as sp

DATA = '/mnt/shared-workspace/dmg_spatial'

REN_SAMPLES = {  # sample -> (gsm, group)
    'DMG1': ('GSM5833528', 'DMG'), 'DMG2': ('GSM5833529', 'DMG'), 'DMG3': ('GSM5833530', 'DMG'),
    'DMG4': ('GSM5833531', 'DMG'), 'DMG5': ('GSM5833532', 'DMG'), 'GBM1': ('GSM5833533', 'GBM'),
    'GBM2': ('GSM5833534', 'GBM'), 'GBM3': ('GSM5833535', 'GBM'), 'GBM4': ('GSM5833536', 'GBM'),
    'GBM5_1': ('GSM5833537', 'GBM'), 'GBM5_2': ('GSM5833538', 'GBM')}

KORD_SAMPLES = {  # region -> (gsm, file_stem)
    'Midbrain': ('GSM8608415', 'VLP73_D1'),
    'PonsA': ('GSM8608416', '_3_VLP70_20221108_D2'),
    'PonsB': ('GSM8608417', '_4_VLP65_20221019_D'),
    'PonsC': ('GSM8608418', 'VLP70_20221108_D1'),
    'PonsD': ('GSM8608419', 'VLP70_20221108_A2'),
    'MedullaA': ('GSM8608420', 'B1_VLP49'),
    'MedullaB': ('GSM8608421', 'C1_VLP45'),
    'MedullaC': ('GSM8608422', 'VLP74_D1'),
    'SpinalCordA': ('GSM8608423', 'VLP74_A1'),
    'SpinalCordB': ('GSM8608424', 'VLP73_A1')}

HUCH_SAMPLES = {'pDMG1': ('GSM8294434', 'pDMG_Sample-1'),
                'pDMG2': ('GSM8294435', 'pDMG_Sample-2'),
                'pDMG3': ('GSM8294436', 'pDMG_Sample-3')}

XENIUM_SAMPLES = {'H3K27M_P6': ('GSM8389375', 'H3K27M_P6_Region_1'),
                  'PFA_P3': ('GSM8389376', 'PFA_P3_Region_3'),
                  'STZR_P2': ('GSM8389377', 'ST-Z_R_P2_Region_1'),
                  'H3WT_hemi_P2': ('GSM8389378', 'H3WT_hemi_P2_Region_1')}


def _finalize(a, sample, cohort, group, positions_path=None, positions_header=0):
    """Attach spatial coordinates and metadata; keep in-tissue spots."""
    a.var_names_make_unique()
    if positions_path is not None:
        pos = pd.read_csv(positions_path, header=positions_header)
        if 'barcode' not in pos.columns:  # header=None case
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
    a.obs['cohort_spatial'] = cohort
    a.obs['group'] = group
    a.X = sp.csr_matrix(a.X)
    return a


def load_ren(sample):
    d = f'{DATA}/ren/extracted/{sample}_spaceranger_out'
    a = sc.read_10x_h5(f'{d}/filtered_feature_bc_matrix.h5')
    return _finalize(a, sample, 'Ren', REN_SAMPLES[sample][1],
                     f'{d}/spatial/tissue_positions_list.csv', positions_header=None)


def load_kordowski(region):
    gsm, stem = KORD_SAMPLES[region]
    d = f'{DATA}/kordowski'
    a = sc.read_10x_h5(f'{d}/{gsm}_{stem}_filtered_feature_bc_matrix.h5')
    return _finalize(a, region, 'Kordowski', 'DMG',
                     f'{d}/{gsm}_{stem}_tissue_positions_list.csv.gz', positions_header=None)


def load_huchede(sample):
    gsm, stem = HUCH_SAMPLES[sample]
    d = f'{DATA}/huchede'
    a = sc.read_mtx(f'{d}/{gsm}_{stem}_matrix.mtx.gz').T
    feats = pd.read_csv(f'{d}/{gsm}_{stem}_features.tsv.gz', sep='\t', header=None)
    a.var_names = feats[1].values if feats.shape[1] > 1 else feats[0].values
    a.obs_names = pd.read_csv(f'{d}/{gsm}_{stem}_barcodes.tsv.gz', sep='\t', header=None)[0].values
    return _finalize(a, sample, 'Huchede', 'DMG',
                     f'{d}/{gsm}_{stem}_tissue_positions.csv.gz', positions_header=0)


def load_all_visium(groups=('DMG', 'GBM', 'peritumor')):
    """Load all Visium sections. Returns dict name -> AnnData."""
    out = {}
    for s, (_, g) in REN_SAMPLES.items():
        if g in groups or (g == 'GBM' and 'GBM' in groups):
            out[f'Ren_{s}'] = load_ren(s)
    for r in KORD_SAMPLES:
        out[f'Kord_{r}'] = load_kordowski(r)
    for s in HUCH_SAMPLES:
        out[f'Huch_{s}'] = load_huchede(s)
    return out


def load_section(name):
    """Load a single section by its full name (e.g. 'Ren_DMG1', 'Kord_PonsA', 'Huch_pDMG1')."""
    if name.startswith('Ren_'):
        return load_ren(name[4:])
    if name.startswith('Kord_'):
        return load_kordowski(name[5:])
    if name.startswith('Huch_'):
        return load_huchede(name[5:])
    if name.startswith('Collot_'):
        from collot_loader import load_collot
        return load_collot(name[7:])
    raise ValueError(name)


if __name__ == '__main__':
    secs = load_all_visium()
    rows = []
    for name, a in secs.items():
        rows.append({'section': name, 'spots': a.shape[0], 'genes': a.shape[1],
                     'median_umi': float(np.median(a.X.sum(1).A1)),
                     'group': a.obs['group'].iloc[0]})
    print(pd.DataFrame(rows).to_string(index=False))

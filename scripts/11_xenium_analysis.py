"""Xenium H3K27M analysis: marker-score state assignment + neighborhood colocalization.

Panel coverage (verified): OPC-like 4/4, AC-like 2/4, OC-like 2/4, Cycling 3/3,
Myeloid 2/4; weak: MES-like 1/4, T cell 1/3, Endothelial 1/3, Neuron 1/2.
"""
import sys, os, json
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import scipy.sparse as sp

DATA = '/mnt/shared-workspace/dmg_spatial'
OUT = '/workspace/dmg_spatial_code/xenium'
os.makedirs(OUT, exist_ok=True)

XEN = {'H3K27M_P6': ('GSM8389375', 'H3K27M_P6_Region_1'),
       'PFA_P3': ('GSM8389376', 'PFA_P3_Region_3'),
       'STZR_P2': ('GSM8389377', 'ST-Z_R_P2_Region_1'),
       'H3WT_hemi_P2': ('GSM8389378', 'H3WT_hemi_P2_Region_1')}

STATE_MARKERS = {'OPC-like': ['PDGFRA', 'OLIG1', 'OLIG2', 'SOX10'],
                 'AC-like': ['GFAP', 'AQP4'],
                 'OC-like': ['MOG', 'MAG'],
                 'MES-like': ['NAMPT'],
                 'Cycling': ['MKI67', 'TOP2A', 'CENPF'],
                 'Myeloid': ['AIF1', 'CD14'],
                 'T cell': ['TRAC'],
                 'Endothelial': ['PECAM1'],
                 'Neuron': ['RBFOX3']}


def load_xenium(sample, with_coords=True):
    gsm, stem = XEN[sample]
    d = f'{DATA}/xenium'
    a = sc.read_mtx(f'{d}/{gsm}_{stem}-matrix.mtx.gz').T
    feats = pd.read_csv(f'{d}/{gsm}_{stem}-features.tsv.gz', sep='\t', header=None)
    a.var_names = feats[1].values
    a.obs_names = pd.read_csv(f'{d}/{gsm}_{stem}-barcodes.tsv.gz', sep='\t', header=None)[0].values
    a.var_names_make_unique()
    a.X = sp.csr_matrix(a.X)
    a.obs['sample'] = sample
    if with_coords:
        import subprocess
        pq = f'/workspace/xenium_tmp/{stem}-transcripts.parquet'
        if not os.path.exists(pq):
            os.makedirs('/workspace/xenium_tmp', exist_ok=True)
            subprocess.run(['bash', '-c',
                            f'gunzip -c {d}/{gsm}_{stem}-transcripts.parquet.gz > {pq}'], check=True)
        t = pd.read_parquet(pq, columns=['cell_id', 'x_location', 'y_location'])
        t = t[t['cell_id'] != 'UNASSIGNED']
        cent = t.groupby('cell_id')[['x_location', 'y_location']].median()
        a.obs = a.obs.join(cent, how='left')
        a.obsm['spatial'] = a.obs[['x_location', 'y_location']].values
        print('cells with coords:', a.obs['x_location'].notna().sum(), '/', a.shape[0])
    return a


def assign_states(a):
    """Score each cell per state with panel markers; assign argmax with confidence margin."""
    sc.pp.normalize_total(a, target_sum=100)  # small target for targeted panel
    sc.pp.log1p(a)
    scores = {}
    for state, gs in STATE_MARKERS.items():
        gs = [g for g in gs if g in a.var_names]
        if gs:
            scores[state] = a[:, gs].X.toarray().mean(1) if sp.issparse(a.X) else a[:, gs].X.mean(1)
    S = pd.DataFrame(scores, index=a.obs_names)
    # z-score each state's score across cells for comparability
    Sz = (S - S.mean()) / S.std()
    a.obs['xenium_state'] = Sz.idxmax(axis=1)
    order = np.argsort(Sz.values, axis=1)
    top = Sz.values[np.arange(len(Sz)), order[:, -1]]
    second = Sz.values[np.arange(len(Sz)), order[:, -2]] if Sz.shape[1] > 1 else 0
    a.obs['xenium_margin'] = top - second
    return a, Sz


if __name__ == '__main__':
    sample = sys.argv[1] if len(sys.argv) > 1 else 'H3K27M_P6'
    a = load_xenium(sample)
    print(f'{sample}: {a.shape[0]} cells x {a.shape[1]} genes')
    a, Sz = assign_states(a)
    print(a.obs['xenium_state'].value_counts())
    print('mean margin:', a.obs['xenium_margin'].mean().round(3))
    a.write(f'{OUT}/{sample}_labeled.h5ad', compression='gzip')
    Sz.to_csv(f'{OUT}/{sample}_state_scores.csv')
    print('saved')

"""cell2location mapping for DMG spatial project — CPU-optimized, shared-output, resumable.

Reference training is complete (ref_signatures.csv on shared volume).
Outputs: /mnt/shared-workspace/c2l_out/<section>/{abundances.csv, qc.json, vis_mapped.h5ad}
"""
import sys, os, time, json, shutil
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')

sys.path.insert(0, '/workspace/dmg_spatial_code')
from visium_loader import load_section

OUT_SHARED = '/mnt/shared-workspace/c2l_out'
OUT_LOCAL = '/workspace/c2l_out'
SIG_PATH = '/mnt/shared-workspace/c2l_out/ref_signatures.csv'
os.makedirs(OUT_LOCAL, exist_ok=True)

CANONICAL = ['PDGFRA','OLIG1','OLIG2','SOX10','GFAP','AQP4','SLC1A3','ALDH1L1','MBP','PLP1','MOG','MAG',
             'CD44','CHI3L1','VIM','NAMPT','MKI67','TOP2A','CENPF','AIF1','CD14','CSF1R','C1QA',
             'CD3D','CD3E','TRAC','VWF','CLDN5','PECAM1','RBFOX3','SYT1']


def discriminative_genes(sig, n_top=2500):
    lv = np.log1p(sig.values)
    vscore = lv.var(axis=1)
    order = np.argsort(-vscore)
    top = list(sig.index[order[:n_top]])
    for g in CANONICAL:
        if g in sig.index and g not in top:
            top.append(g)
    return top


def map_section(name, n_cells=8, detection_alpha=20, max_epochs=2500, n_top_genes=2500):
    import torch
    torch.set_num_threads(16)
    from cell2location.models import Cell2location
    t0 = time.time()
    sig = pd.read_csv(SIG_PATH, index_col=0)
    genes = discriminative_genes(sig, n_top_genes)
    sig = sig.loc[genes]

    vis = load_section(name)
    shared = vis.var_names.intersection(sig.index)
    vis = vis[:, shared].copy()
    sig = sig.loc[shared]
    print(f'[map {name}] spots={vis.shape[0]} genes={len(shared)}', flush=True)

    Cell2location.setup_anndata(vis)
    mod = Cell2location(vis, cell_state_df=sig,
                        N_cells_per_location=n_cells, detection_alpha=detection_alpha)
    mod.train(max_epochs=max_epochs, accelerator='cpu',
              early_stopping=True, early_stopping_monitor='elbo_train',
              early_stopping_patience=30, early_stopping_min_delta=0.01)
    mod.export_posterior(vis, sample_kwargs={'num_samples': 500, 'batch_size': 2500})
    abun = vis.obsm['q05_cell_abundance_w_sf'].copy()
    abun.columns = [c.replace('q05cell_abundance_w_sf_', '') for c in abun.columns]

    sdir_l = f'{OUT_LOCAL}/{name}'; sdir_s = f'{OUT_SHARED}/{name}'
    os.makedirs(sdir_l, exist_ok=True); os.makedirs(sdir_s, exist_ok=True)
    abun.to_csv(f'{sdir_s}/abundances.csv')
    vis.write(f'{sdir_l}/vis_mapped.h5ad', compression='gzip')
    shutil.copy(f'{sdir_l}/vis_mapped.h5ad', f'{sdir_s}/vis_mapped.h5ad')
    hist = mod.history['elbo_train']
    qc = {'section': name, 'spots': int(vis.shape[0]), 'genes': int(len(shared)),
          'minutes': round((time.time()-t0)/60, 1), 'epochs_run': int(len(hist)),
          'elbo_final': float(hist.iloc[-1, 0])}
    with open(f'{sdir_s}/qc.json', 'w') as fh:
        json.dump(qc, fh, indent=2)
    print(f'[map {name}] DONE in {qc["minutes"]} min, {qc["epochs_run"]} epochs', flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'map':
        map_section(sys.argv[2])
    elif sys.argv[1] == 'map_list':
        for name in sys.argv[2:]:
            if os.path.exists(f'{OUT_SHARED}/{name}/abundances.csv'):
                print(f'[skip {name}] done', flush=True)
                continue
            try:
                map_section(name)
            except Exception as e:
                print(f'[map {name}] FAIL: {e}', flush=True)
        print('LIST COMPLETE', flush=True)

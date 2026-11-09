"""Per-section Visium QC + spatial domains via the spatial-transcriptomics skill workflow."""
import sys, os, json, time
import numpy as np, pandas as pd

sys.path.insert(0, '/mnt/skills/system/spatial-transcriptomics/scripts')
sys.path.insert(0, '/workspace/dmg_spatial_code')
from spatial_workflow import run_spatial_analysis
from visium_loader import load_all_visium

OUT = '/workspace/dmg_spatial_code/qc'
os.makedirs(OUT, exist_ok=True)

sections = sys.argv[1:] if len(sys.argv) > 1 else None
if sections:
    from visium_loader import load_section
    secs = {s: load_section(s) for s in sections}
else:
    secs = load_all_visium()
summary = []
for name, a in secs.items():
    if sections and name not in sections:
        continue
    t0 = time.time()
    sdir = f'{OUT}/{name}'
    os.makedirs(sdir, exist_ok=True)
    try:
        a2 = run_spatial_analysis(a, min_genes=200, min_cells=10, max_pct_mito=20,
                                  n_top_genes=2000, n_pcs=30, resolution=0.8,
                                  n_neighbors=15, svgs_n_perms=100, svgs_n_jobs=8,
                                  output_dir=sdir)
        a2.write(f'{sdir}/adata_domains.h5ad', compression='gzip')
        n_svg = int((a2.uns.get('svg_results', pd.DataFrame()).get('fdr', pd.Series([1])) < 0.05).sum()) \
            if 'svg_results' in a2.uns else None
        summary.append({'section': name, 'spots_in': a.shape[0], 'spots_qc': a2.shape[0],
                        'n_leiden': a2.obs['leiden'].nunique() if 'leiden' in a2.obs else None,
                        'n_svg_fdr05': n_svg, 'minutes': round((time.time()-t0)/60, 1), 'status': 'ok'})
        print(f'[qc {name}] ok in {summary[-1]["minutes"]} min', flush=True)
    except Exception as e:
        summary.append({'section': name, 'status': f'FAIL: {e}'})
        print(f'[qc {name}] FAIL: {e}', flush=True)

pd.DataFrame(summary).to_csv(f'{OUT}/qc_summary.csv', index=False)
print('QC SUMMARY')
print(pd.DataFrame(summary).to_string(index=False))

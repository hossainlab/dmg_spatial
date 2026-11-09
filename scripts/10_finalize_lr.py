"""Final LR tables: (1) same-population headline set; (2) interface crosstalk set.

Interface crosstalk evidentiary standard (tiered, per approved plan):
  - recurrent in Kordowski (>=5/10 sections, sig + spatial_z>2), AND
  - significant (magnitude_rank<=0.05) in >=2 Ren DMG sections OR >=1 Huchede section.
"""
import pandas as pd, numpy as np

LR = '/mnt/shared-workspace/lr_out'
OUT = '/mnt/shared-workspace/validation'

REN = ['Ren_DMG1','Ren_DMG2','Ren_DMG3','Ren_DMG4','Ren_DMG5']
KORD = ['Kord_Midbrain','Kord_PonsA','Kord_PonsB','Kord_PonsC','Kord_PonsD',
        'Kord_MedullaA','Kord_MedullaB','Kord_MedullaC','Kord_SpinalCordA','Kord_SpinalCordB']
HUCH = ['Huch_pDMG1','Huch_pDMG2','Huch_pDMG3']
GBM = ['Ren_GBM1','Ren_GBM2','Ren_GBM3','Ren_GBM4','Ren_GBM5_1','Ren_GBM5_2']

def load(secs):
    frames = []
    for s in secs:
        d = pd.read_csv(f'{LR}/{s}_lr_pairs.csv')
        d['section'] = s
        frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    d['interaction'] = d['ligand_complex'] + '|' + d['receptor_complex'] + '|' + d['source'] + '->' + d['target']
    d['cross'] = d['source'] != d['target']
    return d

kord = load(KORD); ren = load(REN); huch = load(HUCH); gbm = load(GBM)

ks = kord[(kord['magnitude_rank'] <= 0.05) & (kord['cross']) & (kord['spatial_z'] > 2)]
krec = ks.groupby('interaction')['section'].nunique()
krec = krec[krec >= 5]
print('Kordowski recurrent interface pairs:', len(krec))

ren_sig = ren[(ren['magnitude_rank'] <= 0.05) & (ren['cross'])]
huch_sig = huch[(huch['magnitude_rank'] <= 0.05) & (huch['cross'])]
gbm_sig = gbm[(gbm['magnitude_rank'] <= 0.05) & (gbm['cross'])]
ren_n = ren_sig.groupby('interaction')['section'].nunique()
huch_n = huch_sig.groupby('interaction')['section'].nunique()
gbm_n = gbm_sig.groupby('interaction')['section'].nunique()

rows = []
for ix, ksec in krec.items():
    rn, hn = int(ren_n.get(ix, 0)), int(huch_n.get(ix, 0))
    if rn >= 2 or hn >= 1:
        sub = ks[ks['interaction'] == ix]
        rows.append({'interaction': ix, 'kord_sections': int(ksec),
                     'ren_dmg_sections_sig': rn, 'huchede_sections_sig': hn,
                     'gbm_sections_sig': int(gbm_n.get(ix, 0)),
                     'mean_lr_means': round(float(sub['lr_means'].mean()), 3),
                     'mean_spatial_z': round(float(sub['spatial_z'].mean()), 2)})
iface = pd.DataFrame(rows).sort_values(['kord_sections','ren_dmg_sections_sig'], ascending=False)
iface.to_csv(f'{OUT}/lr_interface_crosstalk.csv', index=False)
print('Interface crosstalk pairs (Kordowski-recurrent + patient-cohort support):', len(iface))
print(iface.head(25).to_string(index=False))
print('\nGBM support distribution (DMG-selectivity check):')
print(iface['gbm_sections_sig'].value_counts().sort_index())
print('\nDONE')

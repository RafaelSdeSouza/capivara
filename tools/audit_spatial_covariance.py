"""Four-band native covariance diagnostics, not an unvalidated full-spectrum fix."""
from pathlib import Path
import json,csv
import numpy as np
from astropy.io import fits
from scipy.sparse import csr_matrix,diags
P=Path(__file__).resolve().parents[3];O=P/'results/science_baseline_freeze_v3';rows=[]
for path in sorted((O/'physical').glob('*_spaxels.npz')):
 id=path.name.split('_')[0];a=dict(np.load(path));m=a['cluster_map'];ny,nx=m.shape
 with fits.open(P/'sandra-46-trial_bars'/f'manga-{id}-MEGACUBE.fits',memmap=True) as hd:
  for ext in ['GCORREL','RCORREL','ICORREL','ZCORREL']:
   h=hd[ext];t=h.data;bb=float(h.header['BBWAVE']);declared=int(h.header['BBINDEX']);k=int(np.argmin(abs(a['wavelength']-bb)))
   i=t['INDXI_C1']+nx*t['INDXI_C2'];j=t['INDXJ_C1']+nx*t['INDXJ_C2'];rho=np.asarray(t['RHOIJ'],dtype=float)
   assert np.all(i>=0)&np.all(j>=0)&np.all(i<nx*ny)&np.all(j<nx*ny)
   upper=csr_matrix((rho,(i,j)),shape=(nx*ny,nx*ny))
   # The stored triangles contain one unordered pair. Verify before reflecting.
   off=upper-diags(upper.diagonal());assert (off.multiply(off.T)).nnz==0
   corr=upper+upper.T-diags(upper.diagonal())
   for region in np.unique(a['region']):
    ix=np.flatnonzero(a['region']==region);good=np.isfinite(a['variance'][ix,k])&(a['variance'][ix,k]>0)&(a['mask'][ix,k]==0)
    use=ix[good];pixels=a['x'][use]+nx*a['y'][use];std=np.sqrt(a['variance'][use,k]);sub=corr[pixels,:][:,pixels]
    diag=float(std@std);cov=float(std@(sub@std))
    rows.append(dict(galaxy_id=id,region_id=int(region),native_extension=ext,bbwavelength=bb,header_bbindex=declared,
     wave_at_header_bbindex=float(a['wavelength'][declared]),nearest_native_index=k,nearest_native_wavelength=float(a['wavelength'][k]),
     n_region_spaxels=len(ix),n_valid_subset=len(use),variance_diagonal=diag,variance_with_band_spatial_correlation=cov,
     noise_ratio=np.sqrt(cov/diag) if diag>0 else np.nan,
     scope='native MASK==0 subset at BBWAVE; diagnostic only; no spectral covariance or regional fit correction'))
with (O/'variance/native_spatial_correlation_bands.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
print('Band diagnostics',len(rows),'median noise ratio',np.nanmedian([r['noise_ratio'] for r in rows]))
print('Maximum BBINDEX wavelength mismatch',max(abs(r['wave_at_header_bbindex']-r['bbwavelength']) for r in rows))

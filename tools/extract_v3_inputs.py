"""Read-only V2 partitions and native FITS; V3 response-design evidence only."""
from pathlib import Path
import json,hashlib,csv
import numpy as np
from astropy.io import fits
P=Path(__file__).resolve().parents[3];O=P/'results/science_baseline_freeze_v3'
ids=['11004-12701','11014-3704','8602-12705','8932-3701','9869-9102']
meta={r['plateifu']:r for r in csv.DictReader(open(P/'results/science_baseline_freeze_v2/validated/comparison_summary.csv'))}
rows=[]
for id in ids:
 path=P/'sandra-46-trial_bars'/f'manga-{id}-MEGACUBE.fits'
 mp=P/'results/science_baseline_freeze_v2/validated/comparisons'/id/'rest_map.csv'
 m=np.genfromtxt(mp,delimiter=',',skip_header=1).T
 ys,xs=np.where(np.isfinite(m));regions=m[ys,xs].astype(int)
 with fits.open(path,memmap=True,uint=False) as h:
  w=np.asarray(h[6].data,float)
  f=np.asarray(h[1].data[:,ys,xs],float).T
  ivar=np.asarray(h[2].data[:,ys,xs],float).T
  mask=np.asarray(h[3].data[:,ys,xs],int).T
  post=np.asarray(h[4].data[:,ys,xs],float).T
  pre=np.asarray(h[5].data[:,ys,xs],float).T
  assert h[4].name in ['DISP','LSFPOST'] and h[5].name in ['PREDISP','LSFPRE']
  prov=dict(drp_version=h[1].header['VERSDRP3'],native_extension={'pre':h[5].name,'post':h[4].name},
   units='Angstrom',quantity='Gaussian sigma_lambda',wavelength_medium='vacuum',wavelength_frame='observed',
   source_path=str(path),source_sha256=hashlib.file_digest(open(path,'rb'),'sha256').hexdigest(),
   native_bunit={'pre':h[5].header.get('BUNIT'),'post':h[4].header.get('BUNIT')},
   units_source='SDSS native DRP data model; absent container BUNIT',
   orientation='spaxel,wavelength; original Astropy [wavelength,y,x]',
   reference='https://www.sdss4.org/dr17/manga/manga-data/working-with-manga-data/')
  v=np.divide(1.,ivar,out=np.full_like(ivar,np.nan),where=np.isfinite(ivar)&(ivar>0))
  # Entire measured native grid is retained. No fit or segmentation is run here.
  np.savez_compressed(O/'physical'/f'{id}_spaxels.npz',wavelength=w,flux=f,variance=v,mask=mask,
    post=post,pre=pre,region=regions,x=xs,y=ys,cluster_map=m)
  prov.update(galaxy_id=id,redshift=float(meta[id]['redshift']),segmentation_sha256=hashlib.file_digest(open(mp,'rb'),'sha256').hexdigest(),
   segmentation_id='V2-rest-fixed18-a243fdfe',n_spaxels=len(xs))
  (O/'physical'/f'{id}_provenance.json').write_text(json.dumps(prov,indent=2));rows.append(prov)
  print(id,len(xs),flush=True)
(O/'manifests/inputs.json').write_text(json.dumps(rows,indent=2))

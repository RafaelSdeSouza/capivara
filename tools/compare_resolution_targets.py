from pathlib import Path
import json,csv,sys
import numpy as np
from scipy.ndimage import maximum_filter1d
P=Path(__file__).resolve().parents[3];O=P/'results/science_baseline_freeze_v3'
sys.path[:0]=[str(P/'worktrees/capivara-science-baseline-v3/inst/python'),str(P/'worktrees/capivaraPPXF-science-baseline-v3/inst/python')]
from capivara_stellar import template_floor
rows=[];targets=[];data=[]
for path in sorted((O/'physical').glob('*_spaxels.npz')):
 id=path.name.split('_')[0];a=np.load(path);p=json.loads((O/'physical'/f'{id}_provenance.json').read_text());w=a['wavelength']
 k=(w>=3650)&(w<=7800);s=a['post'][:,k];s=np.where(np.isfinite(s)&(s>0),s,np.nan)
 galaxy=np.nanmax(s,axis=0);targets.append(galaxy);data.append((id,a,p,k,s,galaxy))
common=np.nanmax(targets,axis=0)
for id,a,p,k,s,galaxy in data:
 w=a['wavelength'][k];floor=template_floor(w,p['redshift'])
 for reg in np.unique(a['region']):
  x=s[a['region']==reg];mx=np.nanmax(x,axis=0);ref=np.nanmedian(x,axis=0)
  robust=np.nanquantile(x,.95,axis=0)
  for strategy,t in [('local_max',mx),('robust95_unconstrained',robust),('upper_envelope',np.maximum(mx,maximum_filter1d(robust,31,mode='nearest'))),('galaxy_common',galaxy),('five_galaxy_common',common),('local_max_template_floor',np.maximum(mx,floor))]:
   rows.append(dict(galaxy_id=id,region_id=int(reg),strategy=strategy,median_sigma_angstrom=float(np.nanmedian(t)),
     median_sigma_kms=float(np.nanmedian(299792.458*t/w)),median_resolution_ratio=float(np.nanmedian(ref/t)),
     median_added_variance=float(np.nanmedian(t*t-ref*ref)),violating_samples=int(np.sum(t<mx-1e-10)),
     missing_input_samples=int(np.sum(~np.isfinite(x)))))
with (O/'resolution/target_strategies.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
print(len(rows),'strategy rows')

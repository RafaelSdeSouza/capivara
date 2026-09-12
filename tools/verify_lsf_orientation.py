import csv,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]).resolve()
out=root/'results/science_baseline_freeze_v2/physical'
inv={r['id']:r for r in json.loads((out.parent/'lsf_header_inventory.json').read_text())}
rows=list(csv.DictReader(open(out/'lsf_orientation_samples.csv')))
maximum=0
for id in dict.fromkeys(r['plateifu'] for r in rows):
 with fits.open(root/inv[id]['path'],memmap=True) as hdus:
  for r in rows:
   if r['plateifu']!=id:continue
   x,y,k=(int(r[a])-1 for a in ('x','y','channel'))
   for key,hdu in [('flux',1),('post',4),('pre',5)]:
    actual=hdus[hdu].data[k,y,x];expected=float(r[key]) if r[key]!='NA' else np.nan
    if key!='flux' and (not np.isfinite(actual) or actual<=0):actual=np.nan
    maximum=max(maximum,abs(actual-expected))
    assert np.isclose(actual,expected,rtol=1e-12,atol=1e-12,equal_nan=True)
   assert hdus[6].data[k]==float(r['wavelength'])
(out/'orientation_verification.json').write_text(json.dumps({'passed':True,'samples':len(rows),'values':4*len(rows),'max_absolute_difference':maximum,'mapping':'R[x,y,k] equals Astropy[k-1,y-1,x-1] for R one-based indices'},indent=2))
fig,axs=plt.subplots(1,2,figsize=(9,3.5),sharey=True,layout='constrained')
for ax,id in zip(axs,['11004-12701','8602-12705']):
 a=np.genfromtxt(out/(id+'_lsf_profile.csv'),delimiter=',',names=True)
 for key,color,ls in [('pre','#2166ac','-'),('post','#b35806','--')]:
  ax.plot(a['wavelength'],a[key],label=key.upper(),color=color,ls=ls,lw=1.2)
 ax.set(xlabel=r'Observed vacuum wavelength [$\mathrm{\AA}$]',title=f"{id}  ({inv[id]['drp']})")
 ax.set_xlim(3622,10353);ax.tick_params(direction='in',top=True,right=True)
axs[0].set_ylabel(r'LSF $\sigma_\lambda$ [$\mathrm{\AA}$]');axs[0].legend(frameon=False)
fig.savefig(out/'lsf_generations.png',dpi=180)
print('PASS',len(rows),'asymmetric samples;',maximum,'max absolute error')

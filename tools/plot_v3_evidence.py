from pathlib import Path
import csv,json,pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
P=Path(__file__).resolve().parents[3];O=P/'results/science_baseline_freeze_v3';F=O/'figures'
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':170})
def read(path):return list(csv.DictReader(open(O/path)))
def nums(rows,key):return np.array([float(x.get(key,'nan') or 'nan') for x in rows])
def write(path,rows):
 with (O/path).open('w') as f:
  wr=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rows))));wr.writeheader();wr.writerows(rows)
# Scientific parameter recovery, rather than a single global RMSE.
r=read('synthetic/recovery.csv');fig,axs=plt.subplots(2,2,figsize=(8,6),layout='constrained')
for ax,key,label in zip(axs.flat,['delta_log_age','delta_metal','delta_velocity','delta_sigma'],[r'$\Delta\langle\log_{10}t/{\rm yr}\rangle_L$',r'$\Delta\langle[\mathrm{M/H}]\rangle_L$',r'$\Delta v_\star$ [km s$^{-1}$]',r'$\Delta\sigma_\star$ [km s$^{-1}$]']):
 for k,(sigma,color) in enumerate([(20,'#9a6a3a'),(80,'#276a96'),(200,'#555555')]):
  for snr in [10,30,60,150,300]:
   a=[x for x in r if float(x['sigma'])==sigma and float(x['snr_in'])==snr];v=nums(a,key);v=v[np.isfinite(v)]
   q=np.quantile(v,[.16,.5,.84]);ax.errorbar(snr+(k-1)*1.8,q[1],yerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt='o',ms=3,color=color,label=f'{sigma} km/s' if snr==10 else None)
 ax.axhline(0,c='k',lw=.6);ax.set(xlabel='Injected regional S/N before smoothing',ylabel=label);ax.legend(frameon=False,fontsize=7)
fig.savefig(F/'synthetic_recovery.png');plt.close(fig)
line=read('resolution/line_recovery.csv');fig,ax=plt.subplots(figsize=(7,3.5),layout='constrained')
for name in dict.fromkeys(x['shape'] for x in line):
 a=[x for x in line if x['shape']==name];ax.scatter(nums(a,'centre'),100*nums(a,'relative_width_error'),s=10,alpha=.6,label=name)
ax.axhline(0,c='k',lw=.6);ax.axhline(2,c='grey',ls=':',lw=.7);ax.axhline(-2,c='grey',ls=':',lw=.7)
ax.set(xlabel='Observed vacuum wavelength [Angstrom]',ylabel='Recovered width error [%]');ax.legend(frameon=False,ncol=2,fontsize=7)
fig.savefig(F/'resolution_recovery.png');plt.close(fig)
# Full five-pilot inventory; all nominal region IDs remain present.
allrows=[];summary=[]
for folder in sorted((O/'sandra_pilot').iterdir()):
 if not (folder/'regions.csv').exists():continue
 rows=read(str(folder.relative_to(O)/'regions.csv'));assert len(rows)==18
 allrows+=rows;id=folder.name;a=np.load(O/'physical'/f'{id}_spaxels.npz');m=a['cluster_map']
 success=[x for x in rows if x['fit_status'].startswith('SUCCESS')]
 summary.append(dict(galaxy_id=id,n_regions=len(rows),n_numerical_success=len(success),n_bound=sum(x['bound_hit']=='True' for x in rows),
   n_insufficient=sum(x['fit_status']=='INSUFFICIENT_WAVELENGTH' for x in rows),n_numerical_failure=sum(x['fit_status']=='NUMERICAL_FAILURE' for x in rows),
   n_population_usable=sum(x['scientifically_usable']=='True' for x in rows),n_sigma_resolved=sum(x['sigma_resolved']=='True' for x in rows),
   median_snr=float(np.nanmedian(nums(rows,'snr'))),median_chi2=float(np.nanmedian(nums(rows,'chi2'))),
   max_native_flux_relative_L1=float(np.nanmax(nums(rows,'native_relative_L1_error')))))
 # Region map and representative spectra: best residual fit plus an incomplete
 # region where present. Chosen by recorded QC, not population appearance.
 ranked=sorted(success,key=lambda x:float(x['chi2']))
 selected=[int(ranked[0]['region_id'])] if ranked else [1]
 other=[int(x['region_id']) for x in rows if not x['fit_status'].startswith('SUCCESS')]
 selected.append(other[0] if other else int(ranked[-1]['region_id']))
 fig=plt.figure(figsize=(11,6),layout='constrained');grid=fig.add_gridspec(2,3,width_ratios=[1,1.25,1.25]);ax=fig.add_subplot(grid[:,0])
 # A fixed spectral partition is not a galaxy-membership mask. Keep rejected
 # inherited regions visible in grey, without giving them galaxy-region colours.
 # Numerical fitting support still does not certify a stellar population fit.
 supported=np.isin(m,[int(x['region_id']) for x in success])
 ax.imshow(np.where(np.isfinite(m)&~supported,1.,np.nan),origin='lower',
   cmap=ListedColormap(['#e3e3e3']),vmin=0,vmax=1,interpolation='nearest')
 ax.imshow(np.where(supported,m,np.nan),origin='lower',cmap='tab20',vmin=.5,vmax=20.5,interpolation='nearest')
 for region in range(1,19):
  y,x=np.where((m==region)&supported)
  if len(x):
   ax.contour(m==region,levels=[.5],colors='white',linewidths=.3)
   centre=np.argmin((x-np.median(x))**2+(y-np.median(y))**2)
   ax.text(x[centre],y[centre],str(region),ha='center',va='center',fontsize=6,color='black')
 ax.set(xlabel='Native x pixel',ylabel='Native y pixel',title=id+'\nRegions with fitting support')
 ax.legend(handles=[Patch(facecolor='#e3e3e3',label='Insufficient wavelength support')],
   loc='upper center',bbox_to_anchor=(.5,-.16),frameon=False,fontsize=7)
 for j,region in enumerate(selected):
  saved=pickle.load(open(folder/f'region_{region:02d}.pkl','rb'));r=saved['regional'];fit=saved['fit'];p=saved['product']
  ax=fig.add_subplot(grid[j,1]);k=(r['native_wavelength']>3800)&(r['native_wavelength']<7500)
  ax.plot(r['native_wavelength'][k],r['native_flux'][k],color='0.65',lw=.65,label='Native sum')
  ax.plot(r['fitting_wavelength'][k],r['fitting_flux'][k],color='#276a96',lw=.65,label='Controlled spectrum')
  ax.set(xlabel='Observed vacuum wavelength [Angstrom]',ylabel='Summed native flux density',title=f"Region {region}; {fit['fit_status']}")
  ax.legend(frameon=False,fontsize=7)
  ax=fig.add_subplot(grid[j,2])
  if p:
   ax.plot(p['wavelength'],p['galaxy'],lw=.6,color='0.5',label='Fitting spectrum')
   ax.plot(p['wavelength'],p['model'],lw=.65,color='#a65728',label='Stellar model')
   ix=p['goodpixels'];ax.plot(p['wavelength'][ix],p['residual'][ix],lw=.5,color='#276a96',label='Residual')
   ax.set_title(f"Reduced chi-square {fit['chi2']:.2f}; bound={fit['bound_hit']}");ax.legend(frameon=False,fontsize=7)
  else:
   ax.text(.5,.5,'No usable fitting support',transform=ax.transAxes,ha='center');ax.set_axis_off()
  ax.set(xlabel='Rest-air wavelength [Angstrom]',ylabel='Normalized flux density')
 fig.savefig(F/f'{id}_spectroscopy.png');plt.close(fig)
write('sandra_pilot/all_regions.csv',allrows);write('sandra_pilot/pilot_summary.csv',summary)
fig,axs=plt.subplots(2,3,figsize=(10,6),layout='constrained')
snr=nums(allrows,'snr');bad=np.array([x['scientifically_usable']!='True' for x in allrows])
for ax,key,lab in zip(axs[0],['log_age_light','metal_light','sigma_star_raw'],['Mean log age [yr]','Mean template [M/H]','Raw stellar dispersion [km/s]']):
 ax.scatter(snr[bad],nums(allrows,key)[bad],s=12,c='0.6',marker='x',label='QC rejected')
 if (~bad).any():ax.scatter(snr[~bad],nums(allrows,key)[~bad],s=14,c='#276a96',label='Conditional')
 ax.set(xlabel='Formal diagonal S/N',ylabel=lab);ax.legend(frameon=False,fontsize=7)
ratio=nums(allrows,'sigma_star_raw')/nums(allrows,'instrumental_sigma_kms');axs[1,0].scatter(snr,ratio,s=12,c='0.5');axs[1,0].axhline(.8,ls=':',c='k');axs[1,0].set(xlabel='Formal diagonal S/N',ylabel=r'Raw $\sigma_\star/\sigma_{\rm inst}$')
x=np.arange(len(summary));axs[1,1].bar(x-.15,[s['n_bound']/18 for s in summary],.3,label='Bound');axs[1,1].bar(x+.15,[(18-s['n_numerical_success'])/18 for s in summary],.3,label='No numerical fit');axs[1,1].set_xticks(x,[s['galaxy_id'] for s in summary],rotation=45,ha='right',fontsize=7);axs[1,1].set_ylabel('Fraction of 18 regions');axs[1,1].legend(frameon=False,fontsize=7)
chi=nums(allrows,'chi2');axs[1,2].hist(chi[np.isfinite(chi)],bins=12,color='0.6');axs[1,2].axvline(2,c='k',ls=':');axs[1,2].set(xlabel='Reduced chi-square',ylabel='Regions')
fig.savefig(F/'pilot_diagnostics.png');plt.close(fig)
print(json.dumps(summary,indent=2))

"""Label-invariant region correspondence and static science figures for V2.
Usage: python tools/compare_frame_products.py VALIDATED_OUTPUT
"""
import os, sys, json, csv
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/capivara-v2-matplotlib')
import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]); figdir=root/'figures';figdir.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
summary=[]
for d in sorted((root/'comparisons').iterdir()):
    if not (d/'comparison.rds').exists(): continue
    maps={k:np.genfromtxt(d/f'{k}_map.csv',delimiter=',',skip_header=1) for k in ['observed','rest']}
    sums={k:np.genfromtxt(d/f'{k}_regional_sums.csv',delimiter=',',skip_header=1) for k in maps}
    reg={k:np.genfromtxt(d/f'{k}_regions.csv',delimiter=',',names=True) for k in maps}
    t=np.genfromtxt(d/'contingency.csv',delimiter=',',skip_header=1)
    a,b=linear_sum_assignment(t,maximize=True)
    # One-to-one maximum-overlap assignment is a comparison convention only;
    # retain the full contingency matrix to expose splitting and merging.
    rows=[]
    for i,j in zip(a,b):
        fa,fb=sums['observed'][:,i+1],sums['rest'][:,j+1]
        good=np.isfinite(fa)&np.isfinite(fb)
        sa,sb=np.sum(np.abs(fa[good])),np.sum(np.abs(fb[good]))
        l1=np.sum(np.abs(fb[good]-fa[good]))/sa if sa else np.nan
        shape=np.sum(np.abs(fb[good]/sb-fa[good]/sa)) if sa and sb else np.nan
        ra,rb=reg['observed'][i],reg['rest'][j]
        rows.append(dict(observed_region=i+1,rest_region=j+1,overlap_spaxels=int(t[i,j]),
            observed_size=int(ra['n_spaxels']),rest_size=int(rb['n_spaxels']),
            size_delta=int(rb['n_spaxels']-ra['n_spaxels']),region_iou=t[i,j]/(ra['n_spaxels']+rb['n_spaxels']-t[i,j]),
            observed_bar_fraction=ra['bar_fraction'],rest_bar_fraction=rb['bar_fraction'],
            bar_fraction_delta=rb['bar_fraction']-ra['bar_fraction'],
            observed_components=int(ra['n_components']),rest_components=int(rb['n_components']),
            summed_spectrum_relative_L1=l1,unit_L1_spectrum_difference=shape,
            compared_native_channels=int(good.sum())))
    with (d/'matched_region_comparison.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    def values(k): return np.array([r[k] for r in rows])
    summary.append(dict(plateifu=d.name,matched_spaxel_fraction=float(t[a,b].sum()/t.sum()),
        median_abs_size_change=float(np.median(np.abs(values('size_delta')))),
        median_abs_bar_fraction_change=float(np.median(np.abs(values('bar_fraction_delta')))),
        max_abs_bar_fraction_change=float(np.max(np.abs(values('bar_fraction_delta')))),
        median_summed_spectrum_relative_L1=float(np.nanmedian(values('summed_spectrum_relative_L1'))),
        max_summed_spectrum_relative_L1=float(np.nanmax(values('summed_spectrum_relative_L1'))),
        median_unit_L1_spectrum_difference=float(np.nanmedian(values('unit_L1_spectrum_difference')))))
    matched=np.full_like(maps['rest'],np.nan)
    for i,j in zip(a,b):matched[maps['rest']==j+1]=i+1
    fig,axs=plt.subplots(1,2,figsize=(7,3.4),layout='constrained')
    for ax,m,label in zip(axs,[maps['observed'],matched],['4800–7400 Å observed','4800–7400 Å rest']):
        ax.imshow(m,origin='lower',interpolation='nearest',cmap='tab20',vmin=.5,vmax=20.5)
        ax.set(xlabel='Column',ylabel='Row',title=label)
    fig.savefig(figdir/f'{d.name}_partitions.png');plt.close(fig)
    # Full summed spectra for the four largest historical regions. The paired
    # candidate is chosen solely by spatial overlap, never by spectral similarity.
    largest=np.argsort([r['observed_size'] for r in rows])[-4:][::-1]
    fig,axs=plt.subplots(2,2,figsize=(8,5),layout='constrained')
    for ax,k in zip(axs.flat,largest):
        r=rows[k];i,j=r['observed_region'],r['rest_region'];wave=sums['observed'][:,0]
        ax.plot(wave,sums['observed'][:,i],lw=.7,color='#336699',label=f'Observed region {i}')
        ax.plot(wave,sums['rest'][:,j],lw=.7,alpha=.85,color='#aa5500',label=f'Rest region {j}')
        ax.set(xlabel='Native wavelength [Å]',ylabel='Summed native $F_\lambda$')
        ax.legend(fontsize=7,frameon=False)
    fig.savefig(figdir/f'{d.name}_summed_spectra.png');plt.close(fig)
with (root/'regional_comparison_summary.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
(root/'regional_metrics_definitions.json').write_text(json.dumps({
  'correspondence':'One-to-one maximum native-spaxel overlap assignment (SciPy linear_sum_assignment); contingency retained',
  'summed_spectrum_relative_L1':'sum(abs(S_rest-S_observed))/sum(abs(S_observed)) on common finite native wavelengths',
  'unit_L1_spectrum_difference':'L1 distance after each matched sum spectrum is divided by its own L1 norm; dimensionless',
  'spectra':'Full retained native FLUX sums; no new mask or IVAR cut; not fitted or interpolated',
  'partitions':'Native row/column grid, origin lower. Colours matched by maximum spatial overlap, not an assertion of region identity.',
  'spectral_panels':'Four largest historical regions and their spatial-overlap-matched candidates; full native summed flux in original units, no per-region normalisation.'},indent=2))
print(json.dumps(summary,indent=2))

"""Five frozen partitions only; refuses to run before synthetic gate completion."""
from pathlib import Path
import sys,json,csv,pickle,time
import numpy as np
P=Path(__file__).resolve().parents[3];O=P/'results/science_baseline_freeze_v3'
sys.path[:0]=[str(P/'worktrees/capivara-science-baseline-v3/inst/python'),str(P/'worktrees/capivaraPPXF-science-baseline-v3/inst/python')]
from capivara_resolution import prepare_region
from capivara_stellar import StellarLibrary,fit_region,template_floor
G=json.loads((O/'synthetic/gates.json').read_text());assert G['allow_pilot'],'Synthetic gate is closed'
LIB=StellarLibrary(P/'outputs/capivara_ppxf_manga_smoke/sps_models/spectra_emiles_9.0.npz')
ids=['11004-12701','11014-3704','8602-12705','8932-3701','9869-9102']
if len(sys.argv)>1:
 assert sys.argv[1] in ids;ids=[sys.argv[1]]
for id in ids:
 a=np.load(O/'physical'/f'{id}_spaxels.npz');prov=json.loads((O/'physical'/f'{id}_provenance.json').read_text());out=O/'sandra_pilot'/id;out.mkdir(exist_ok=True)
 rows=[];fluxcheck=[];w=a['wavelength'];z=prov['redshift'];floor=template_floor(w,z)
 for region in np.unique(a['region']):
  choose=a['region']==region;f=a['flux'][choose];v=a['variance'][choose];post=a['post'][choose];pre=a['pre'][choose]
  validlsf=np.isfinite(post)&(post>0);safe=np.where(validlsf,post,1.)
  target=np.maximum(safe.max(0),floor)
  saved=out/f'region_{int(region):02d}.pkl'
  if saved.exists():
   r=pickle.load(open(saved,'rb'))['regional']
  else:
   r=prepare_region(w,f,v,pre,post,a['mask'][choose]==0,id,int(region),prov['segmentation_id'],z,prov,
     strategy='common',target=target,retain_operators=False)
  r['configuration']['strategy']='regional maximum with explicit template-resolution floor'
  # Compare native full-grid sums against immutable V2 measured sums.
  v2=np.genfromtxt(P/'results/science_baseline_freeze_v2/validated/comparisons'/id/'rest_regional_sums.csv',delimiter=',',skip_header=1)
  old=v2[:,int(region)];both=np.isfinite(old)&np.isfinite(r['native_flux'])
  relative=np.sum(abs(old[both]-r['native_flux'][both]))/np.sum(abs(old[both]))
  assert relative<1e-12
  result,product,config=fit_region(r,LIB)
  resolved=(np.isfinite(result['sigma_star_raw']) and result['snr']>=30 and result['sigma_star_raw']/result.get('instrumental_sigma_kms',np.inf)>=.8 and 65<=result['sigma_star_raw']<=220 and not result['bound_hit'])
  eligible=(result['fit_status']=='SUCCESS_INTERIOR' and result['snr']>=30 and result['chi2']<=2 and not result['population_edge'])
  if eligible:
   result,product,config=fit_region(r,LIB,bootstrap=20,seed=20260911+int(region))
  result['scientifically_usable']=bool(eligible and np.isfinite(result['log_age_uncertainty']))
  result['sigma_resolved']=bool(resolved and result['scientifically_usable'])
  result['sigma_star']=result['sigma_star_raw'] if result['sigma_resolved'] else np.nan
  reason=[]
  if result['snr']<30:reason.append('BELOW_VALIDATED_SNR')
  if result['chi2']>2:reason.append('EXCESS_RESIDUAL_OR_UNDERSTATED_VARIANCE')
  if result['population_edge']:reason.append('TEMPLATE_GRID_EDGE')
  if not resolved:reason.append('DISPERSION_OUTSIDE_VALIDATED_DOMAIN')
  if result['bound_hit']:reason.append('PARAMETER_BOUND')
  if not result['fit_status'].startswith('SUCCESS'):reason.append(result['reason'])
  result['qc_restrictions']=';'.join(reason)
  if result['scientifically_usable']:result['reason']='CONDITIONAL_SAME_FAMILY_DOMAIN_SPATIAL_COVARIANCE_UNKNOWN'
  result.update(native_relative_L1_error=float(relative),fitting_masked_fraction=r['masked_fraction'],
    selected_lsf='POST',target_sigma_median=float(np.median(target)),drp=prov['drp_version'],
    n_input_spaxels=int(choose.sum()),valid_fitting_channels=int(r['valid'].sum()),
    min_fitting_native=float(w[r['valid']].min()) if np.any(r['valid']) else np.nan,
    max_fitting_native=float(w[r['valid']].max()) if np.any(r['valid']) else np.nan)
  pickle.dump(dict(regional=r,fit=result,product=product,configuration=config),open(out/f'region_{int(region):02d}.pkl','wb'),protocol=4)
  tab={k:v for k,v in result.items() if k!='bootstrap_draws'};rows.append(tab)
  with (out/'regions.csv').open('w') as h:
   wr=csv.DictWriter(h,fieldnames=sorted(set().union(*(r.keys() for r in rows))));wr.writeheader();wr.writerows(rows)
  print(id,int(region),result['fit_status'],'snr',round(result['snr'],1),'chi2',round(result['chi2'],2),'usable',result['scientifically_usable'],flush=True)
 print('PILOT COMPLETE',id,len(rows),flush=True)

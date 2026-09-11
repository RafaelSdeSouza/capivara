from pathlib import Path
import csv,json
import numpy as np
P=Path(__file__).resolve().parents[3];O=P/'results/science_baseline_freeze_v3'
def read(p):return list(csv.DictReader(open(O/p)))
def write(p,r):
 with (O/p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(r[0]));w.writeheader();w.writerows(r)
r=read('synthetic/recovery.csv');assert len(r)==810,'Incomplete recovery grid'
summary=[]
for key in ['snr_in','sigma','age','metal','drp']:
 for val in sorted(set(x[key] for x in r)):
  a=[x for x in r if x[key]==val];d=dict(stratum=key,value=val,n=len(a),failed=sum(not x['fit_status'].startswith('SUCCESS') for x in a),bound_fraction=np.mean([x['bound_hit']=='True' for x in a]))
  for par in ['log_age','metal','velocity','sigma']:
   v=np.array([float(x['delta_'+par]) for x in a]);v=v[np.isfinite(v)]
   d[par+'_bias']=float(np.median(v));d[par+'_scatter']=float(1.4826*np.median(abs(v-np.median(v))));d[par+'_p90_absolute_error']=float(np.quantile(abs(v),.9))
  summary.append(d)
write('synthetic/stratified_recovery.csv',summary)
# Conditional domain: require every age/metal/DRP marginal within the high-SNR,
# resolved-injection subset. Marginal tests do not certify every interpolation.
domain=[x for x in r if float(x['snr_in'])>=30 and float(x['sigma'])>=80]
gates=[]
for key in ['age','metal','drp','sigma','snr_in']:
 for val in sorted(set(x[key] for x in domain)):
  a=[x for x in domain if x[key]==val];d=dict(stratum=key,value=val,n=len(a),pass_gate=True)
  for par,tol,scattertol in [('log_age',.15,.25),('metal',.15,.2),('velocity',10,20),('sigma',15,25)]:
   v=np.array([float(x['delta_'+par]) for x in a]);v=v[np.isfinite(v)]
   bias=float(np.median(v));scatter=float(1.4826*np.median(abs(v-bias)))
   d[par+'_bias']=bias;d[par+'_scatter']=scatter;d['pass_gate'] &= abs(bias)<=tol and scatter<=scattertol
  d['pass_gate'] &= np.mean([x['fit_status'].startswith('SUCCESS') for x in a])>=.95
  gates.append(d)
write('synthetic/domain_gates.csv',gates)
resolution=read('resolution/line_recovery.csv');variance=read('variance/monte_carlo.csv');cov=read('variance/parameter_covariance_comparison.csv')
cp=[]
for key in set((x['galaxy_id'],x['sigma'],x['rep']) for x in cov):
 a=[x for x in cov if (x['galaxy_id'],x['sigma'],x['rep'])==key]
 assert len(a)==2
 a.sort(key=lambda x:x['covariance_treatment']);d=dict(galaxy_id=key[0],sigma=key[1],rep=key[2])
 for par in ['log_age_light','metal_light','v_star','sigma_star_raw']:
  d['full_minus_diagonal_'+par]=float(a[1][par])-float(a[0][par])
 cp.append(d)
write('variance/paired_parameter_differences.csv',cp)
checks=dict(resolution=max(abs(float(x['relative_width_error'])) for x in resolution)<.02,
 variance=max(abs(float(x['variance_relative_median'])) for x in variance)<.02,
 covariance_all_success=all(x['fit_status'].startswith('SUCCESS') for x in cov),
 recovery_domain=all(x['pass_gate'] for x in gates),
 covariance_parameter_agreement=all(np.nanmedian([abs(x['full_minus_diagonal_'+p]) for x in cp])<tol for p,tol in [('log_age_light',.05),('metal_light',.05),('v_star',5),('sigma_star_raw',5)]))
verdict=dict(checks=checks,allow_pilot=all(checks.values()),
 regional_resolution_contract='CONDITIONAL',variance_contract='DIAGONAL_APPROXIMATION_VALIDATED' if checks['variance'] and checks['covariance_all_success'] else 'FAILED',
 stellar_population_recovery='CONDITIONAL_SNR30_AGE0.1to10_METALminus0.71to0.22_SIGMA80to200' if checks['recovery_domain'] else 'FAILED',
 native_error_assumption='independent pixels and spaxels; retained induced spectral covariance',
 population_basis='36 E-MILES templates; conditional marginal recovery, no arbitrary SFH identification',
 sigma_reporting_rule='SNR_fitting>=30 and sigma_raw/instrumental_sigma>=0.8, no bound hit; domain limited to injected 80--200 km/s',
 uncertainty='bootstrap SD conditional on fitted spectrum; no model or spatial covariance uncertainty')
(O/'synthetic/gates.json').write_text(json.dumps(verdict,indent=2));print(json.dumps(verdict,indent=2))

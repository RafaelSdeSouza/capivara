"""Frozen same-family grid, mixtures, scalar comparison and mismatch controls."""
import os,sys,json,csv,pickle,itertools
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import CubicSpline
P=Path(__file__).resolve().parents[3];O=P/'results/science_baseline_freeze_v3'
sys.path[:0]=[str(P/'worktrees/capivara-science-baseline-v3/inst/python'),str(P/'worktrees/capivaraPPXF-science-baseline-v3/inst/python')]
from capivara_resolution import gaussian_operator,prepare_region
from capivara_stellar import StellarLibrary,fit_region,air_to_vacuum,template_floor,C
LIB=StellarLibrary(P/'outputs/capivara_ppxf_manga_smoke/sps_models/spectra_emiles_9.0.npz')
IDS=['11004-12701','8602-12705']
W=np.arange(3650.,7801.)
# Fixed before parameter recovery. Every row and failed fit is retained.
BASE_GRID=list(itertools.product(IDS,[.1,1.,10.],[-.71,0.,.22],[20.,80.,200.],[10.,30.,60.],range(3)))
GRID=BASE_GRID+list(itertools.product(IDS,[.1,1.,10.],[-.71,0.,.22],[20.,80.,200.],[150.,300.],range(3)))
CONFIG=dict(seed=20260911,objects=IDS,ages_gyr=[.1,1.,10.],metals=[-.71,0.,.22],sigma_kms=[20,80,200],snr=[10,30,60,150,300],seeds_per_cell=3,
 velocity_kms=[-120,0,120],n_grid=len(GRID),basis_ages=LIB.age.tolist(),basis_metals=LIB.metal.tolist(),template_sha256=LIB.sha256,
 criteria={'median_bias_logage_dex':.15,'median_bias_metal_dex':.15,'median_bias_velocity_kms':10,'median_bias_sigma_kms':15,
 'robust_scatter_logage_dex':.25,'robust_scatter_metal_dex':.2,'max_bound_fraction':.1},
 limitation='Same-family conditional grid, three realizations/cell; not an independent SSP-model validation')
if __name__=='__main__' and (len(sys.argv)<2 or sys.argv[1] in ['grid','smoke','extended']):
 (O/'synthetic/configuration.json').write_text(json.dumps(CONFIG,indent=2))
profiles={id:np.genfromtxt(P/'results/science_baseline_freeze_v2/physical'/f'{id}_lsf_profile.csv',delimiter=',',names=True) for id in IDS}
provs={id:json.loads((O/'physical'/f'{id}_provenance.json').read_text()) for id in IDS}

def truth_flux(age,metal,vel,sigma,z,lsf):
 # Forward model starts on an independent fine log grid, Gaussian LOSVD then
 # native wavelength-dependent instrumental broadening in quadrature with
 # the known empirical template resolution. No pPXF call generates truth.
 j=np.argmin(abs(np.log(LIB.age/age))+abs(LIB.metal-metal))
 lam=LIB.wave;dv=5.;ln=np.arange(np.log(lam[0]),np.log(lam[-1]),dv/C)
 high=CubicSpline(lam,LIB.flux[:,j])(np.exp(ln));high=gaussian_filter1d(high,sigma/dv,mode='nearest')
 shifted=CubicSpline(ln,high)(np.log(lam)-vel/C)
 air=lam;vac=air_to_vacuum(air);obs=vac*(1+z)
 # Native instrument broadening before sampling; interpolation to a finer
 # 0.1 A observed grid controls the forward-model quadrature error.
 fine=np.arange(W[0]-30,W[-1]+30,.1)
 g=CubicSpline(obs,shifted)(fine)
 inst=np.interp(fine,W,lsf);templ=np.interp(fine,obs,LIB.sigma*(1+z)*vac/air)
 additional=np.sqrt(np.maximum(0,inst**2-templ**2))
 # scipy variable output-centred fine-grid quadrature, independent of the
 # discrete second-moment-calibrated production operator.
 half=int(np.ceil(6*additional.max()/.1));out=np.zeros_like(g);den=np.zeros_like(g)
 for k in range(-half,half+1):
  weight=np.exp(-.5*(k*.1/np.maximum(additional,.001))**2)
  vals=np.interp(fine+k*.1,fine,g);out+=weight*vals;den+=weight
 out/=den
 return CubicSpline(fine,out)(W)

cache={};rows=[]
def run_case(id,age,metal,sigma,snr,rep,mixture=False,approx='homogenized',mismatch=False,full=False):
 z=provs[id]['redshift'];p=profiles[id];base=np.interp(W,p['wavelength'],p['post'])
 # Keep templates resolvable at every wavelength; this forward grid does not
 # simulate deconvolution of unresolved features of the empirical library.
 base=np.maximum(base,template_floor(W,z)*1.03)
 lsf=np.array([base,base*1.12,base*1.25]);pre=np.sqrt(lsf**2-.12)
 velocity=[-120.,0.,120.][rep%3]
 key=(id,age,metal,sigma,rep%3,mixture)
 if key not in cache:
  fs=[]
  for j,s in enumerate(lsf):
   aa=age if not mixture else [age,1.,10.][j]
   zz=metal if not mixture else [metal,0.,-.4][j]
   fs.append(truth_flux(aa,zz,velocity,sigma,z,s))
  cache[key]=np.array(fs)
 noiseless=cache[key].copy()
 if mismatch:
  # Localized abundance-response perturbation, deliberately outside the basis.
  for centre in [5175.,5270.]:
   observed=float(air_to_vacuum(centre)*(1+z));noiseless*=1-.04*np.exp(-.5*((W-observed)/7)**2)
 rng=np.random.default_rng(20260911+GRID.index((id,age,metal,sigma,snr,rep)) if (id,age,metal,sigma,snr,rep) in GRID else 20262000+rep)
 # SNR refers to the regional sum before smoothing at unit mean continuum.
 err=np.median(noiseless,axis=1)*np.sqrt(3)/snr
 variance=np.broadcast_to(err[:,None]**2,noiseless.shape).copy()
 flux=noiseless+rng.normal(size=noiseless.shape)*np.sqrt(variance)
 target=lsf.max(0)
 r=prepare_region(W,flux,variance,pre,lsf,np.ones_like(flux,bool),id,1,'synthetic',z,provs[id],retain_operators=False)
 if approx!='homogenized':
  from scipy.sparse import diags
  r['fitting_flux']=r['native_flux'];r['fitting_variance']=r['native_variance'];r['covariance']=diags(r['native_variance']).tocsr()
  r['target_sigma']=np.full(W.size,np.median(lsf)) if approx=='median_scalar' else np.sqrt(np.mean(lsf**2,axis=0))
  r['valid']=np.ones(W.size,bool);r['valid'][:30]=False;r['valid'][-30:]=False
 result,product,conf=fit_region(r,LIB,full_covariance=full)
 truthage=np.log10(age*1e9) if not mixture else np.mean(np.log10(np.array([age,1.,10.])*1e9))
 truthmetal=metal if not mixture else np.mean([metal,0.,-.4])
 result.update(age=age,metal=metal,sigma=sigma,snr_in=snr,rep=rep,drp=provs[id]['drp_version'],approximation=approx,
  mixture=mixture,mismatch=mismatch,true_velocity=velocity,true_log_age=truthage,true_metal=truthmetal,
  delta_log_age=result['log_age_light']-truthage,delta_metal=result['metal_light']-truthmetal,
  delta_velocity=result['v_star']-velocity,delta_sigma=result['sigma_star_raw']-sigma)
 return result,product,conf

def write(path,rr):
 with (O/path).open('w') as f:
  wr=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rr))));wr.writeheader();wr.writerows(rr)
if __name__=='__main__':
 mode=sys.argv[1] if len(sys.argv)>1 else 'grid'
 if mode=='smoke':
  result,product,conf=run_case(IDS[0],1.,0.,80.,60.,1)
  print(result,flush=True);pickle.dump((result,product,conf),open(O/'synthetic/smoke.pkl','wb'));sys.exit(0)
 if mode in ['grid','extended']:
  if mode=='extended':rows=list(csv.DictReader(open(O/'synthetic/recovery.csv')));assert len(rows)==486
  for n,case in enumerate(GRID if mode=='grid' else GRID[486:]):
   r,prod,c=run_case(*case);rows.append(r)
   if n%20==0:print(n,len(GRID),r['fit_status'],r['delta_log_age'],r['delta_sigma'],flush=True);write('synthetic/recovery.csv',rows)
  write('synthetic/recovery.csv',rows)
 elif mode=='comparisons':
  for id,age,sigma,rep,mix,approx in itertools.product(IDS,[.1,10.],[20.,80.,200.],range(3),[False,True],['homogenized','median_scalar','rms_variable']):
   r,prod,c=run_case(id,age,0.,sigma,60.,rep,mixture=mix,approx=approx);rows.append(r)
  write('synthetic/regional_mixtures.csv',rows)
 elif mode=='mismatch':
  for id,age,sigma,rep,mis in itertools.product(IDS,[1.,10.],[80.,200.],range(3),[False,True]):
   r,prod,c=run_case(id,age,0.,sigma,60.,rep,mismatch=mis);rows.append(r)
  write('template_mismatch/recovery.csv',rows)
 elif mode=='covariance':
  for id,sigma,rep,full in itertools.product(IDS,[20.,80.,200.],range(3),[False,True]):
   r,prod,c=run_case(id,1.,0.,sigma,30.,rep,full=full);rows.append(r)
   write('variance/parameter_covariance_comparison.csv',rows)
 print('FINISHED',mode,len(rows),flush=True)

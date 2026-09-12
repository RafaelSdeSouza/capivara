"""Independent response and covariance experiments; never writes V1/V2."""
import os,sys,json,csv,hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'inst/python'))
from capivara_resolution import gaussian_operator,prepare_region,resolution_target
P=Path(__file__).resolve().parents[3];O=P/'results/science_baseline_freeze_v3'
def savecsv(path,rows):
 with (O/path).open('w') as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
def gauss(x,a,mu,s):return a*np.exp(-.5*((x-mu)/s)**2)
def width(w,f,centre):
 k=abs(w-centre)<35
 pars=curve_fit(gauss,w[k],f[k],p0=[f.max(),centre,2],maxfev=5000)[0]
 return abs(pars[2])
# High-resolution quadrature of intrinsic Gaussian through a variable response.
# This does not use the production smoothing operator or its moment calibration.
def forward(w,centre,intrinsic,lsf):
 x=np.linspace(centre-7*intrinsic,centre+7*intrinsic,1401)
 weights=np.exp(-.5*((x-centre)/intrinsic)**2);weights/=weights.sum()
 s=np.interp(x,w,lsf);k=abs(w-centre)<50
 f=np.zeros(w.size)
 f[k]=(np.exp(-.5*((w[k,None]-x)/s)**2)/(np.sqrt(2*np.pi)*s) @ weights)
 return f
w=np.arange(3622.,10001.);rows=[]
shapes={'constant':np.full(w.size,1.3),'smooth':1.1+.00013*(w-w[0]),
 'sharp':1.4+.16*np.tanh((w-6000)/10)}
for id in ['11004-12701','8602-12705']:
 p=np.genfromtxt(P/'results/science_baseline_freeze_v2/physical'/f'{id}_lsf_profile.csv',delimiter=',',names=True)
 shapes[id]=p['post'][p['wavelength']<=10000]
for name,base in shapes.items():
 for scale in [.8,1.]:
  s=base*scale;t=base*1.05
  a=gaussian_operator(w,np.sqrt(t*t-s*s))
  for centre in [4000.,5000.,5990.,6000.,6010.,7000.]:
   for intrinsic in [.15,1.,5.]:
    f=forward(w,centre,intrinsic,s);goal=forward(w,centre,intrinsic,t);out=a@f
    sig=width(w,out,centre);truth=width(w,goal,centre)
    rows.append(dict(shape=name,scale=scale,centre=centre,intrinsic_sigma=intrinsic,
      sigma_recovered=sig,sigma_target=truth,relative_width_error=sig/truth-1,
      recovered_sigma_kms=sig/centre*299792.458,target_sigma_kms=truth/centre*299792.458,
      flux_error=out.sum()/f.sum()-1,profile_relative_L2=np.linalg.norm(out-goal)/np.linalg.norm(goal)))
savecsv('resolution/line_recovery.csv',rows)
# Variance and induced covariance, including subpixel kernels and no-op.
rng=np.random.default_rng(20260911);small=np.arange(4800.,4928.);vr=[]
for sigma in [0.,.2,.5,1.,2.]:
 a=gaussian_operator(small,np.full(small.size,sigma));v=np.linspace(.4,1.4,small.size)
 cov=(a@ __import__('scipy').sparse.diags(v)@a.T).toarray()
 noise=rng.normal(size=(20000,small.size))*np.sqrt(v);draw=noise@a.T
 empirical=np.cov(draw,rowvar=False)
 vr.append(dict(convolution_sigma=sigma,variance_relative_median=np.median(np.diag(empirical)/np.diag(cov)-1),
   covariance_relative_frobenius=np.linalg.norm(empirical-cov)/np.linalg.norm(cov),
   lag1_correlation=np.median(np.diag(cov,1)/np.sqrt(np.diag(cov)[:-1]*np.diag(cov)[1:])),
   continuum_error_underestimate=np.sqrt(cov.sum()/np.trace(cov))))
 if sigma==0:assert np.array_equal(a.toarray(),np.eye(small.size))
savecsv('variance/monte_carlo.csv',vr)
meta=dict(wavelength_grid=[float(w[0]),float(w[-1]),float(w[1]-w[0])],seed=20260911,
 acceptance={'width_relative_tolerance':.02,'flux_relative_tolerance':1e-12,'profile_L2_tolerance':.02,
 'variance_median_tolerance':.02,'native_LSF_local_variation_mask':.05},
 model='Independent fine-grid Gaussian response quadrature; production discrete variable operator applied only afterward',
 scope='Gaussian native POST response premise; no gas intrinsic correction')
(O/'resolution/experiment.json').write_text(json.dumps(meta,indent=2))
print('Worst width',max(abs(r['relative_width_error']) for r in rows))
print(vr)

"""Regional response experiment: native sums and a separate fitting product.

Gaussian POST LSFs and independent native pixel errors are assumptions, not
measurements established by this module. Variable kernels preserve column flux;
no single Gaussian response is attached to the native sum.
"""
import hashlib
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.optimize import brentq
from scipy.ndimage import maximum_filter1d


def gaussian_operator(wave, sigma):
    """Column-normalized variable Gaussian, discrete second-moment calibrated.

Uniform native grids only. Sigma=0 is exactly identity. Moment calibration
avoids the almost-identity error of sampling a subpixel Gaussian. It does not
assert that a three-pixel kernel is itself Gaussian. Edges are flagged by the
caller. The 6-sigma truncation is part of the stored contract.
    """
    w=np.asarray(wave,float); s=np.asarray(sigma,float)
    if w.ndim!=1 or w.size<3 or s.shape!=w.shape or not np.all(np.isfinite(w)):
        raise ValueError('Invalid wavelength/sigma dimensions')
    step=np.diff(w)[0]
    if step<=0 or not np.allclose(np.diff(w),step,rtol=1e-8,atol=1e-8):
        raise ValueError('Homogenization requires a uniform native grid')
    if np.any(~np.isfinite(s)) or np.any(s<0):raise ValueError('Missing/negative convolution sigma')
    pix=s/step; radius=max(1,int(np.ceil(6*np.max(pix))))
    offsets=np.arange(-radius,radius+1)
    # A tabulated inverse is deterministic and accurate to 1e-5 pixel variance.
    # Zero is kept exact; the nonzero branch includes kernels narrower than a pixel.
    qgrid=np.geomspace(.045,max(1.,float(np.max(pix))*1.1),4096)
    weights=np.exp(-.5*(offsets[:,None]/qgrid[None,:])**2)
    moments=(offsets[:,None]**2*weights).sum(0)/weights.sum(0)
    q=np.interp(pix**2,moments,qgrid)
    rows=[];cols=[];data=[]
    for off in offsets:
        j=np.arange(w.size);i=j+off;ok=(i>=0)&(i<w.size)
        v=np.exp(-.5*(off/q)**2)
        v[pix==0]=float(off==0)
        rows.extend(i[ok]);cols.extend(j[ok]);data.extend(v[ok])
    a=csr_matrix((data,(rows,cols)),shape=(w.size,w.size))
    a=a @ diags(1/np.asarray(a.sum(axis=0)).ravel())
    a.eliminate_zeros()
    return a.tocsr()


def resolution_target(lsf, strategy='local_max', common=None):
    s=np.asarray(lsf,float)
    if s.ndim!=2 or np.any(~np.isfinite(s)) or np.any(s<=0):
        raise ValueError('Target design requires complete positive LSF coverage')
    exact=s.max(axis=0)
    if strategy=='local_max':return exact
    if strategy=='upper_envelope':
        # A robust percentile alone violates the contract. Constrain it by max.
        return np.maximum(exact,maximum_filter1d(np.quantile(s,.95,axis=0),size=31,mode='nearest'))
    if strategy=='common':
        t=np.asarray(common,float)
        if t.shape!=exact.shape or np.any(~np.isfinite(t)) or np.any(t<exact-1e-12):
            raise ValueError('Common target is narrower than a contributor')
        return t
    raise ValueError('Unknown target strategy')


def prepare_region(wavelength, flux, variance, lsf_pre, lsf_post, mask,
                   galaxy_id, segment_id, segmentation_id, redshift,
                   lsf_provenance, strategy='local_max', target=None,
                   retain_operators=False):
    """A strict fixed-spaxel regional sum; incomplete fitting channels are masked.

The native finite-FLUX sum is unchanged by IVAR/MASK cuts. Native variance is
missing unless every finite-flux contributor has positive finite variance.
Fitting channels require all contributors and all kernel support to be valid.
    """
    w=np.asarray(wavelength,float);f=np.atleast_2d(np.asarray(flux,float));v=np.asarray(variance,float)
    pre=np.asarray(lsf_pre,float);post=np.asarray(lsf_post,float);m=np.asarray(mask,bool)
    if any(x.shape!=f.shape for x in [v,pre,post,m]) or f.shape[1]!=w.size:
        raise ValueError('Regional arrays must all be (spaxel,wavelength)')
    if not np.isfinite(redshift) or redshift<=-1:raise ValueError('Invalid systemic redshift')
    if lsf_provenance.get('wavelength_medium')!='vacuum' or lsf_provenance.get('wavelength_frame')!='observed':
        raise ValueError('This stellar preparation contract requires observed vacuum wavelengths')
    if lsf_provenance.get('quantity')!='Gaussian sigma_lambda' or lsf_provenance.get('units')!='Angstrom':
        raise ValueError('Native LSF units/quantity are required')
    finite=np.isfinite(f);native=np.where(finite,f,0).sum(0); counts=finite.sum(0)
    native[counts==0]=np.nan
    vok=np.isfinite(v)&(v>0)
    native_var=np.where(finite&vok,v,0).sum(0)
    native_var[np.any(finite&~vok,axis=0)|(counts==0)]=np.nan
    goodlsf=np.isfinite(post)&(post>0)
    full=goodlsf.all(0)
    # No interpolation over missing LSF. A placeholder only constructs the
    # operator; all affected outputs are explicitly invalidated below.
    clean=np.where(goodlsf,post,1.)
    t=resolution_target(clean,strategy,target)
    conv=np.sqrt(np.maximum(0,t[None,:]**2-clean**2))
    fit=np.zeros(w.size);cov=csr_matrix((w.size,w.size));fit_counts=np.zeros(w.size,int);ops=[]
    for j in range(f.shape[0]):
        a=gaussian_operator(w,conv[j]);support=a.copy();support.data[:]=1
        valid=finite[j]&vok[j]&m[j]&goodlsf[j]&full
        complete=np.asarray(support @ (~valid).astype(float)).ravel()==0
        # Conservative halo: no renormalized edge profile is accepted.
        radius=max(1,int(np.ceil(6*np.max(conv[j])/np.diff(w)[0])))
        complete[:radius]=False;complete[-radius:]=False
        fit+=np.asarray(a @ np.where(valid,f[j],0)).ravel()
        cov+=a @ diags(np.where(valid,v[j],0)) @ a.T
        fit_counts+=complete
        if retain_operators:ops.append(a)
    valid=fit_counts==f.shape[0]
    # Abrupt native-response changes are outside the locally Gaussian premise.
    radius=max(1,int(np.ceil(6*np.max(t)/np.diff(w)[0])))
    variation=np.zeros(w.size)
    for s in clean:
        hi=maximum_filter1d(s,size=2*radius+1,mode='nearest')
        lo=-maximum_filter1d(-s,size=2*radius+1,mode='nearest')
        variation=np.maximum(variation,(hi-lo)/s)
    valid &= (variation<=.05) & (np.min(clean,axis=0)/np.diff(w)[0]>=.9)
    fitting_var=cov.diagonal();valid &= np.isfinite(fitting_var)&(fitting_var>0)
    fitting_flux=fit.copy();fitting_flux[~valid]=np.nan;fitting_var[~valid]=np.nan
    return dict(schema='capivara_regional_spectrum_v3_candidate',galaxy_id=str(galaxy_id),
      segment_id=int(segment_id),segmentation_id=str(segmentation_id),n_spaxels=f.shape[0],
      native_wavelength=w,native_flux=native,native_variance=native_var,native_counts=counts,
      fitting_wavelength=w.copy(),fitting_flux=fitting_flux,fitting_variance=fitting_var,
      fitting_counts=fit_counts,valid=valid,redshift=float(redshift),
      input_wavelength_frame='observed',fitting_wavelength_frame='observed',wavelength_medium='vacuum',
      lsf_pre=pre,lsf_post=post,lsf_provenance=lsf_provenance,selected_lsf='post',
      target_sigma=t,convolution_sigma=conv,combination='unweighted sum after per-spaxel smoothing',
      native_combination='finite native FLUX sum; no Gaussian effective LSF',
      masked_fraction=float(1-valid.mean()),resolution_variation=variation,
      flux_conservation_status='NATIVE_EXACT_FINITE_SUM',resolution_status='CONDITIONAL_GAUSSIAN_POST',
      variance_status='SPECTRAL_COVARIANCE_PROPAGATED_SPATIAL_INDEPENDENCE_ASSUMED',
      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      qc='PREPARED' if valid.sum()>=50 else 'INSUFFICIENT_WAVELENGTH',
      covariance=cov,operators=ops,configuration=dict(strategy=strategy,kernel='column Gaussian discrete second-moment calibrated',
      truncation_sigma=6,maximum_lsf_fractional_variation=.05,minimum_native_lsf_sigma_pixels=.9,model_pixel_sampling='point_sampled',
      template_pixel_integration=False,gas_fitting=False,missing_policy='strict fixed contributors and kernel halo'))

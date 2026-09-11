import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'inst/python'))
from capivara_resolution import gaussian_operator,resolution_target,prepare_region
class ResponseTests(unittest.TestCase):
 def setUp(self):
  self.w=np.arange(4800.,5101.);self.f=np.ones((2,len(self.w)));self.v=np.ones_like(self.f)
  self.s=np.array([np.full(len(self.w),1.1),np.full(len(self.w),1.5)])
  self.prov=dict(quantity='Gaussian sigma_lambda',units='Angstrom',wavelength_medium='vacuum',wavelength_frame='observed',drp_version='v3_1_1',native_extension={'pre':'LSFPRE','post':'LSFPOST'})
 def prepare(self,**kw):
  args=dict(wavelength=self.w,flux=self.f,variance=self.v,lsf_pre=self.s*.95,lsf_post=self.s,mask=np.ones_like(self.f,bool),galaxy_id='test',segment_id=4,segmentation_id='immutable',redshift=.03,lsf_provenance=self.prov,retain_operators=True)
  args.update(kw);return prepare_region(**args)
 def test_noop(self):
  a=gaussian_operator(self.w,np.zeros(len(self.w)));np.testing.assert_array_equal(a.toarray(),np.eye(len(self.w)))
 def test_native_and_covariance(self):
  r=self.prepare();np.testing.assert_array_equal(r['native_flux'],self.f.sum(0));np.testing.assert_array_equal(r['native_variance'],self.v.sum(0))
  expected=sum((a@a.T for a in r['operators']));np.testing.assert_allclose(r['covariance'].toarray(),expected.toarray())
  self.assertTrue(np.all(r['target_sigma']>=self.s.max(0)));self.assertFalse(r['valid'][0]);self.assertEqual(r['segment_id'],4)
 def test_missing_not_imputed(self):
  s=self.s.copy();s[0,100]=np.nan;r=self.prepare(lsf_post=s)
  self.assertFalse(r['valid'][100]);self.assertTrue(np.isnan(r['fitting_flux'][100]));self.assertEqual(r['native_flux'][100],2)
 def test_invalid_variance_and_mask(self):
  v=self.v.copy();v[0,80]=0;r=self.prepare(variance=v)
  self.assertTrue(np.isnan(r['native_variance'][80]));self.assertFalse(r['valid'][80])
  m=np.ones_like(self.f,bool);m[:,100]=False;r=self.prepare(mask=m);self.assertFalse(r['valid'][100])
 def test_generations(self):
  a=self.prepare();p=self.prov.copy();p.update(drp_version='v2_7_1',native_extension={'pre':'PREDISP','post':'DISP'})
  b=self.prepare(lsf_provenance=p);np.testing.assert_array_equal(a['fitting_flux'],b['fitting_flux']);self.assertEqual(b['lsf_provenance']['native_extension']['post'],'DISP')
 def test_guards(self):
  with self.assertRaisesRegex(ValueError,'narrower'):resolution_target(self.s,'common',np.ones(len(self.w)))
  with self.assertRaises(ValueError):gaussian_operator(self.w,np.full(len(self.w),np.nan))
  with self.assertRaisesRegex(ValueError,'uniform'):gaussian_operator(self.w**2,np.ones(len(self.w)))
  p=self.prov.copy();p['wavelength_medium']='air'
  with self.assertRaisesRegex(ValueError,'vacuum'):self.prepare(lsf_provenance=p)
 def test_second_moment(self):
  for sigma in [.1,.3,1.,2.]:
   a=gaussian_operator(self.w,np.full(len(self.w),sigma));f=a[:,150].toarray().ravel()
   self.assertLess(abs(np.sum(f*(self.w-self.w[150])**2)-sigma**2),2e-5)
 def test_failed_region_retained(self):
  r=self.prepare(mask=np.zeros_like(self.f,bool));self.assertEqual(r['qc'],'INSUFFICIENT_WAVELENGTH');self.assertEqual(r['segment_id'],4)
if __name__=='__main__':unittest.main()

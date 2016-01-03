func circularPupFunction(x,y,obs){
  r = abs(x,y);
  return r<=1.0 & r>obs;
}

func OTF_telescope(void)
/* DOCUMENT otf = OTF_telescope()

   Computes the OTF of the telescope, so that
   > fft(OTF_telescope()).re
   produces a PSF normalized with max(psf)=SR=1.0
     
   SEE ALSO:
 */
{
  N = data.fourier.npix;
  // computation of pupil
  x = data.fourier.ud / (data.tel.diam/2.) * (indgen(N)-(N/2+1));   // x exprime en rayon pupille

  pup =  circularPupFunction(x(,-),x(-,),data.tel.obs);
  
  // factor that will normalize the psf
  // with PSF(0)=1.00 when diffraction-limited
  surface_pup_m2 = data.tel.diam^2*(1-data.tel.obs^2)*pi/4;
  surface_pup_pix = surface_pup_m2 / data.fourier.ud^2;
  factnorm = surface_pup_pix^2;
  // compute FTO of telescope. Computing the psf using
  // just fft(FTO).re produces a psf with max(psf)=SR
  // In fact, FTOtel is normalized so that sum(FTOtel)=1.
  FTOtel = autocorrelation(pup) / factnorm;

  return FTOtel;
}

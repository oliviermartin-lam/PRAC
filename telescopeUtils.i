func circularPupFunction(D,obs,N,pixSize)
/* DOCUMENT

 */
{
  //defining geometry
  x     = (indgen(N) -N/2)(,-:1:N);
  x    *= pixSize/(D/2.);
  y     = transpose(x);
    
  r = abs(x,y);
  return r<=1.0 & r>obs;
}

func polarFromCartesian(x,y)
/* DOCUMENT
 */
{
  rho   = abs(x,y);
  theta = 0*rho;

  if(numberof(x) == 1){
    theta = atan(y/x);
    if(x>0 & y<0){
      theta += 2*pi;
    }else if(x<0){
      theta += pi;
    }else if(x==0 & y>0){
      theta = pi/2.;
    }else if(x==0 & y<0){
      theta = 3*pi/2.;
    }
  }else{
    w0        = where(x>0);
    theta(*)(w0) = atan(y(*)(w0)/x(*)(w0));
    w0        = where(x>0 & y<0);
    theta(*)(w0) = atan(y(*)(w0)/x(*)(w0)) + 2*pi;
    w0        = where(x<0);
    theta(*)(w0) = atan(y(*)(w0)/x(*)(w0)) + pi;
    w0        = where(x==0 & y>0);
    theta(*)(w0) = pi/2.;
    w0        = where(x==0 & y<0);
    theta(*)(w0) = 3*pi/2.;
  }
  return [rho,theta];
}

func OTF_telescope(D,obs,nPixels,pixSize)
/* DOCUMENT otf = OTF_telescope()

   Computes the OTF of the telescope, so that
   > fft(OTF_telescope()).re
   produces a PSF normalized with max(psf)=SR=1.0
     
   SEE ALSO:
 */
{
  // computation of pupil
  pup =  circularPupFunction(D,obs,nPixels,pixSize);
  OTF_tel = fft(abs(fft(pup))^2).re / numberof(pup)^2;
  // Normalize the OTF to get a psf with PSF(0)=1.00 when diffraction-limited
  //OTF_tel /= sum(OTF_tel);
  
  return OTF_tel;
}

func computeOTFstatic(&PSF_stats)
{

  //defining geometry
  N     = tel.nPixels;
  x     = (indgen(N) -N/2)(,-:1:N);
  x    *= tel.pixSize/(tel.diam/2.);
  y     = transpose(x);

  //computes static aberrations in arcsec
  stats  = (*rtc.slopes_res)(slrange(rtc.its),avg);
  
  //computing the Zernike modes in radians
  a  = (*sys.slopesToZernikeMatrix)(,+) * stats(+);
  a *=  pi*tel.diam/(radian2arcsec*cam.lambda);

  //Generating and saving zernike modes
  Zi = readfits("fitsFiles/zernikeModes.fits",err=1);
  if(is_void(Zi)){
    tmp   = polarFromCartesian(x,y);
    rho   = tmp(,,1);
    theta = tmp(,,2);
    Zi = array(0.0,N,N,35);
    for(i=1;i<=35;i++){
      Zi(,,i) = computeZernikePolynomials(rho,theta,i);
    }
    writefits,"fitsFiles/zernikeModes.fits",Zi;
  }

  //loop on Zernike modes
  phi  = array(0.,N,N);
  for(i=3;i<=35;i++){
    phi += a(i) * Zi(,,i) ;
  }
  
  //FTO of the telescope + abstats
  P = circularPupFunction(tel.diam,tel.obs,N,tel.pixSize);
  OTF_stats = autocorrelation(P*exp(1i*phi));

  //normalization
  surface_pup_m2 = tel.diam^2*(1-tel.obs^2)*pi/4;
  surface_pup_pix = surface_pup_m2 / tel.pixSize^2;
  factnorm = surface_pup_pix^2;
  OTF_stats /= factnorm;

  PSF_stats = fft(OTF_stats).re;
  
  return OTF_stats;
}


/*

 _   _  ____ ____   _    
| \ | |/ ___|  _ \ / \   
|  \| | |   | |_) / _ \  
| |\  | |___|  __/ ___ \ 
|_| \_|\____|_| /_/   \_\
                         
*/
func computeNcpaPsd(N,size,pitch,strehl,power,lambda)
/* DOCUMENT

   
 */
{
  //Spatial frequency array
  k      = computeSpatialFreqRad(N, size, kx, ky);
  msk    = defineDmFrequencyArea(k, kx, ky, "square", pitch );
  northo = where(!msk);
    
  //converting the best bench SR into a wavefront variance using the Marechal approx.
  WFE_ncpa     = sqrt(- (lambda*1e9/2/pi)^2 * log(strehl)); // in nm

  //determining the PSD energy in rd^2
  ncpaEnergy   = (2*pi* WFE_ncpa / size / (lambda*1e9))^2; //in rd^2

  
  //The static abberations in the instrument are modelled by a k^-power spectrum
  ncpaSpectrum = 0*k;
  ncpaSpectrum(northo) = (k(northo))^(-power);
  ncpaSpectrum *= ncpaEnergy/sum(ncpaSpectrum);

  return ncpaSpectrum;
}

func calibrateNcpaPsd(&param)
/* DOCUMENT

 */
{

  //Getting the real bench ncpa OTF
  DirNCPA = "2013_09_15_onsky/";
  nPx = cam.nPixelsCropped;
  otfNCPA =  getOTFncpa(nPx,DirNCPA);

   
  //Grabbing the telescope OTF without static terms
  size = 2*tel.diam/nPx;
  otfTel  = OTF_telescope(tel.diam,tel.obs,nPx,size);

  // Managing the fitting
  inputs = array(pointer,4);
  inputs(1) = &otfTel;
  inputs(2) = &tel.fourierPixSize;
  inputs(3) = &dm.pitch;
  inputs(4) = &cam.lambda;

  param = [.9,2];
  error;
  res   = lmfit(modelOtfNcpa,inputs,param,roll(otfNCPA));

  return modelOtfNcpa(inputs,param);
}


func modelOtfNcpa(inputs,a)
/* DOCUMENT

 */
{
  otfTel = *inputs(1);
  N      = dimsof(otfTel)(0);
  size   = *inputs(2);
  pitch  = *inputs(3);
  lambda = *inputs(4);
  
  // Getting the PSD
  strehl = abs(a(1));
  power  = a(2); 
  PSD = computeNcpaPsd(N,size,pitch,strehl,power,lambda);

  // Transforming the PSF into an OTF
  PSD(N/2+1,N/2+1) = 0.0;
  PSD(N/2+1,N/2+1) = -sum(PSD);  
  DPHI = 2*abs(fft(PSD)) *  (size^2);
  
  return otfTel * exp(-0.5*DPHI);
  

}
/*
  
 ____  _____ ____ ___  _   ___     _____  _    _   _ _____ ___ ___  _   _ 
|  _ \| ____/ ___/ _ \| \ | \ \   / / _ \| |  | | | |_   _|_ _/ _ \| \ | |
| | | |  _|| |  | | | |  \| |\ \ / / | | | |  | | | | | |  | | | | |  \| |
| |_| | |__| |__| |_| | |\  | \ V /| |_| | |__| |_| | | |  | | |_| | |\  |
|____/|_____\____\___/|_| \_|  \_/  \___/|_____\___/  |_| |___\___/|_| \_|
                                                                          
*/

func deconvolveFromOTFtelescope(OTF,OTF_tel,method,param)
/*DOCUMENT

 */
{

  OTFold = OTFnew = OTF;
  
  if(method == "raw"){
    
    OTFnew = OTFold/(OTF_tel + param);
       
  }else if(method == "SVD"){
  
    //performing svd of OTF_tel
    l = SVdec(OTF_tel, U, V);
    //OTF_tel = (U*l(-,))(,+)*V(+,)

    //inversion
    l1 = invert_eigen( l, nmf );
    //deconvolution
    w = where(abs(log10(l)(dif)) == max(abs(log10(l)(dif))))(1);
    l1 = 0*l;
    l1(1:w) = 1./l(1:w);
    H = (U*l1(-,))(,+)*V(+,)
    OTFnew = H * OTFold;

  }else if(method =="Bracewell"){
    
    H = 0*OTF_tel;
    w = where(OTF_tel/max(OTF_tel) > param);
    (H(*))(w) = 1./(OTF_tel(*))(w);
    
    OTFnew = H * OTFold;
    
  }else if(method == "Wiener"){

    Q = param;
    K = transpose(OTF_tel)/(OTF_tel*transpose(OTF_tel) + Q^2);
    OTFnew = abs(K*Q)^2 + abs(OTFold - K*OTF_tel)^2;
   
  }
  

 return OTFnew;
}













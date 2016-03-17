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

func circularAveragePsf(im)
/*DOCUMENT
 */
{
  // Geometry
  N   = dimsof(im)(0);
  x   = (indgen(N) -N/2-1)(,-:1:N);
  y   = transpose(x);
  tmp = polarFromCartesian(x,y);
  rho = tmp(,,1);
  t   = tmp(,,2); 
  rmax = max(rho);
  
  // Initialisation
  nt = 32;nr = N/2;//nr = 128;
  imAvg = npx = array(0.,nr,nt);
 
  for(i=1;i<=nt;i++){
    if(i == 1){
      msk_t = t < 2*pi/double(nt);
    }else{
      msk_t = (t >= (i-1)*2*pi/double(nt) ) & (t < i*2*pi/double(nt));
    }
    for(k=1;k<=nr;k++){
      if(k == 1){
        msk_r = rho < rmax/double(nr);
      }else{
        msk_r = rho >= (k-1)*rmax/double(nr) & rho < k*rmax/double(nr);
      }

      mskIm = im * msk_t * msk_r;
      w     = where(mskIm !=0);
      if(is_array(w)){
        imAvg(k,i) = avg(mskIm(w));
        npx(k,i)   = 1.;
      }
    }
  }

  w     = where(npx(,sum) !=0);
  imAvg = imAvg(w,sum)/npx(w,sum);
  
  return imAvg;
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
  OTF_tel /= sum(OTF_tel);
  
  return OTF_tel;
}

func computeOTFstatic(&PSF_stats,&PSF_ncpa_fit,&SR_stats,&SR_ncpa,&sigNcpa)
{

  //defining geometry
  N     = tel.nPixels;
  x     = (indgen(N) -N/2-1)(,-:1:N);
  x    *= tel.pixSize/(tel.diam/2.);
  y     = transpose(x);


  //computes static aberrations in arcsec
  res = (*rtc.slopes_res)(slrange(rtc.its),);
  stats  = res(,avg);
  
  //computing the Zernike modes in radians
  a  = (*sys.slopesToZernikeMatrix)(,+) * stats(+);
  a  = grow(0.,a); // adding piston;
  a *=  pi*tel.diam/(radian2arcsec*cam.lambda);
  
  //Grabbing ncpa calib.
  nmodesmax = 150;
  ancpa  = readfits("fitsFiles/ncpaCalibCoefs_"+var2str(nmodesmax)+".fits"); //in rd
  nmodes = numberof(ancpa);
  mod    = readfits("fitsFiles/ncpaCalibModes_"+var2str(nmodesmax)+".fits"); //in rd
  npx    = dimsof(mod)(2); 
    
  //Generating and saving zernike modes
  Zi = readfits("fitsFiles/zernikeModes.fits",err=1);
    
  if(is_void(Zi)){
    tmp   = polarFromCartesian(x,y);
    rho   = tmp(,,1);
    theta = tmp(,,2);
    Zi = array(0.0,N,N,nmodes);
    for(i=1;i<=nmodes;i++){
      Zi(,,i) = computeZernikePolynomials(rho,theta,i); // i=1 -> piston
    }
    writefits,"fitsFiles/zernikeModes.fits",Zi;
  }
    

  //loop on Zernike modes for TT
  phi_tt  = array(0.,N,N);
  for(i=2;i<=3;i++){
    phi_tt += a(i) * Zi(,,i);
  }
  
  //loop on Zernike modes for CPA
  phi  = array(0.,N,N);
  for(i=4;i<=36;i++){
    phi += a(i) * Zi(,,i);
  }


  //loop on Zernike modes for NCPA
  phi_n  = array(0.,npx,npx);
  for(i=4;i<=nmodes;i++){
    phi_n += ancpa(i) * mod(,,i);
  }

  
  if(npx !=N){
    //Interpolate NCPA phase
    phi_ncpa = 0*phi;
    x0 = (indgen(npx) - npx/2+1)(,-:1:npx) *tel.pixSize * N/(npx*tel.diam/2.);
    y0 = transpose(x0);
    phi_ncpa = interp2(y,x,phi_n,y0,x0);
  }else{
    phi_ncpa = phi_n;
  }
  
  //FTO of the telescope + abstats
  P = circularPupFunction(tel.diam,tel.obs,N,tel.pixSize);
  OTF_stats = autocorrelation(P*exp(1i*(phi + phi_ncpa + phi_tt)));

  //normalization
  OTF_stats /= tel.aeraInPix^2;
  PSF_stats = fft(OTF_stats).re;

    
  //computing statistics
  SR_stats = sum(autocorrelation(P*exp(1i*(phi))))/sum(autocorrelation(P));
  SR_ncpa  = sum(autocorrelation(P*exp(1i*(phi_ncpa))))/sum(autocorrelation(P));
  sigNcpa  = phi_ncpa(where(P))(rms) * 1550./2./pi;

  //Getting the fitting NCPA PSF
  PSF_ncpa_fit = roll(abs(fft(roll(P*exp(1i*phi_ncpa))))^2);
  PSF_ncpa_fit /= max(PSF_ncpa_fit)/SR_ncpa;

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













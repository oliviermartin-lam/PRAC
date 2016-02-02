

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
  OTF_tel = fft(abs(fft(pup))^2).re;
  // Normalize the OTF to get a psf with PSF(0)=1.00 when diffraction-limited
  OTF_tel /= sum(OTF_tel);
  
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



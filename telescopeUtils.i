func residualOTF(method,verb=,disp=)
/* DOCUMENT residualOTF,"instantaneous",verb=,disp=1

 */
{

  //defining geometry
  N     = tel.nPixels;
  x     = (indgen(N) -N/2)(,-:1:N);
  x    *= tel.pixSize/(tel.diam/2.);
  y     = transpose(x);
  

  //computes residues in arcsec
  residue  = (*rtc.slopes_res)(slrange(rtc.its),);
  stats    = residue(,avg);
  dyn      = residue - stats;
  
  if(method == "instantaneous"){
    
    ////////////////////////////////////////////////////////
    // .... Determining the denoised and dealiased residual slopes
    ////////////////////////////////////////////////////////////////

    // .... MMSE reconstrutor
    //empirical covariance matrix of the residue
    C_all = dyn(,+) * dyn(,+)/dimsof(dyn)(0);
    //inversion
    C_all_m = invgen(C_all,0,cond=30.);
    //closed-loop noise matrix
    varNoise  = getNoiseVar(dyn);//in arcsec^2
    Cnn       = 0*C_all;
    takesDiag(Cnn) = varNoise;
    //aliasing matrix
    if(rtc.obsMode == "MOAO"){
      Crr    = *covMatrix.aliasing;
      Crr    = computeCovSlopesError(Crr,*rtc.R);
    } else if(rtc.obsMode == "SCAO"){
      Crr = (*covMatrix.aliasing)(slrange(rtc.its),slrange(rtc.its));
    }
    C_ee = C_all - Cnn  - Crr;
    //MMSE reconstruction
    Rmmse  = C_ee(,+) * C_all_m(+,);
    est    = Rmmse(,+) * dyn(+,); // in arcsec
    
    
    ////////////////////////////////////////////////////
    // .... Average the Zernike modes every exposure time
    //////////////////////////////////////////////////////////////

    zer2rad = pi*tel.diam/(radian2arcsec*cam.lambda); 
    //computing the Zernike modes in radians
    a_dyn  = (*sys.slopesToZernikeMatrix)(,+) * est(+,) *zer2rad;
    a_stat = (*sys.slopesToZernikeMatrix)(,+) * stats(+)*zer2rad;

    // ... Generating and saving zernike modes
    Zi = readfits("fitsFiles/zernikeModes.fits",err=1);
    if(is_void(Zi)){
      Zi = array(0.0,N,N,35);
      for(i=1;i<=35;i++){
        Zi(,,i) = computeZernikePolynomials(rho,theta,i);
      }
      writefits,"fitsFiles/zernikeModes.fits",Zi;
    }
    
    //constants
    t       = int(cam.exposureTime*rtc.Fe);
    nframes = dimsof(dyn)(0);
    nstep   = int(nframes/t);
    nzer = dimsof(a_dyn)(2);
    a_res_avg = a_stat_avg = array(0.,nzer,nstep+1);

    // Loop on exposure time step
    for(i = 1;i<=nstep;i++){
      tk = 1 + (i-1)*t:i*t;//kth exposition
      //average on dynamic modes: the phase is averaged
      a_res_avg(,i)  = (a_dyn(,tk))(,avg);
      //static modes: the slopes are average to get the static phase
      stat_k = (residue(,tk))(,avg);
      a_stat_k = (*sys.slopesToZernikeMatrix)(,+) * stat_k(+);
      a_stat_avg(,i) = a_stat_k*zer2rad;
    }
    
    // .... managing the last frames
    //static
    stat_k = (residue(,nstep*t:))(,avg);//arcsec
    a_stat_k = (*sys.slopesToZernikeMatrix)(,+) * stat_k(+);
    a_stat_avg(,0) = a_stat_k*zer2rad;
    //dynamic
    a_res_avg(,0) =(a_dyn(,t*nstep:))(,avg);

    //Adding static 
    a_res_avg    += a_stat_avg;

    //////////////////////////////
    // .... Determining the PSF
    ////////////////////////////////////////////
    P = circularPupFunction(x,y,tel.obs);
    phi_res_k = PSF_res = 0*P;
    
    for(k=1;k<=nstep+1;k++){
      phi_res_k *=0;
      for(i=1;i<=35;i++){
        phi_res_k += a_res_avg(i,k) * Zi(,,i);
      }
      //Electrical Field
      E = roll(P*exp(1i*phi_res_k));
      //PSF
      PSF_res += abs(fft(E))^2;
      write,format="Job done:%.3g%s\r",100.*k/nstep,"%";
    }

    /////////////////////////////////
    // .... Adding high order modes and normalization
    /////////////////////////////////////////////////////
    
    // .... Adding fitting + ncpa
    OTF_res = fft(PSF_res).re;
    OTF_res *= roll(*otf.fit)*computeFakeOTFncpa(budget.SRncpa);

    //grabbing the PSF
    PSF_res = roll(fft(OTF_res).re);
    PSF_res /= sum(PSF_res);
    PSF_res /= tel.airyPeak;
    
    SR_res = arrondi(100*max(PSF_res),1);

    //cropping the reconstructed PSF
    nm = (tel.nPixels - cam.nPixelsCropped)/2+1;
    np = (tel.nPixels +  cam.nPixelsCropped)/2;
    PSF_res =  PSF_res(nm:np,nm:np);
    
  } else if(method == "Veran"){
    // .... Computing covariance matrices
    stat     = residue(,avg);
    residue -= stat;
    //TS measurements covariance matrix
    Cee  = residue(,+) * residue(,+)/dimsof(residue)(0);
    // noise matrix
    Cnn  = 0*Cee;
    takesDiag(Cnn) = varNoise(residue);
    // aliasing matrix
    Crr = *covMatrix.aliasing;
    Crr = computeCovSlopesError(Crr,*rtc.R);
    //true residual covariance matrix estimation
    Cres = Cee - Cnn + Crr;
    
    Dphi_res = 2*(Cres - max(Cres)); 
    OTF_res  = exp(-0.5*Dphi_res);

  }


  /////////////////////
  // ..... Getting Ensquared Energy on both reconstructed/sky PSF
  //////////////////////////////////////////////////////////
  
  boxsize = EE  = span(1,cam.nPixelsCropped-3.,cam.nPixelsCropped) * cam.pixSize;
  PSF_sky = *psf.sky;
   
  for(i=1; i<=cam.nPixelsCropped; i++){
    EE(i) = getEE( 100*PSF_res/sum(PSF_res), cam.pixSize, boxsize(i));
  }
  
  ///////////////
  // .... Display
  ///////////////////////////////

  if(disp){
    window,6;clr;
    N = cam.nPixelsCropped;
    pli,PSF_res,-N*tel.pixSize,-N*tel.pixSize,N*tel.pixSize,N*tel.pixSize;
    pltitle,"Direct reconstructed PSF with SR = "+var2str(SR_res)+"%";
    xytitles,"Arcsec","Arcsec";
  
    window,7;clr;
    pli,*psf.sky,-N*tel.pixSize,-N*tel.pixSize,N*tel.pixSize,N*tel.pixSize;
    pltitle,"Sky PSF with SR = "+var2str(arrondi(100*max(*psf.sky),1))+"%";
    xytitles,"Arcsec","Arcsec";
  
    winkill,8;window,8,style="aanda.gs",dpi=90;clr;
    plg, *psf.EE_sky, boxsize;
    plg, *psf.EE_res,boxsize,color=[100,100,100];
    plg, EE, boxsize,color=[128,128,128];
    plg, [100,100],[-0.1,max(boxsize)*1.05],type=2,marks=0;
    fcut = radian2arcsec*cam.lambda/tel.pitch;
    plg, [0,100],[fcut,fcut],type=2,marks=0;
    plg, [0,100],[1,1]*1.22*radian2arcsec*cam.lambda/tel.diam,type=2,marks=0;
    xytitles,"Angular separation from center (arcsec)","Ensquared Energy (%)";
    plt,"A: On-sky PSF",1.5,40,tosys=1;
    plt,"B: OTF-based reconstructed PSF",1.5,30,tosys=1;
    plt,"C: Telemetry-based reconstructed PSF",1.5,20,tosys=1;
    plt,"DM cut frequency",fcut*1.05,50,tosys=1;
    plt,"1.22 !l/D",1.22*radian2arcsec*cam.lambda/tel.diam*1.05,10,tosys=1;
    range,0,105;
    limits,-0.1,max(boxsize)*1.05;
  }

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
  x   =  pixSize / (D/2.) * (indgen(nPixels)-(nPixels/2+1));
  pup =  circularPupFunction(x(,-),x(-,),obs);
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
  tmp   = polarFromCartesian(x,y);
  rho   = tmp(,,1);
  theta = tmp(,,2);

  //computes static aberrations in arcsec
  stats  = (*rtc.slopes_res)(slrange(rtc.its),avg);
  
  //computing the Zernike modes in radians
  a  = (*sys.slopesToZernikeMatrix)(,+) * stats(+);
  a *=  pi*tel.diam/(radian2arcsec*cam.lambda);

  //Generating and saving zernike modes
  Zi = readfits("fitsFiles/zernikeModes.fits",err=1);
  if(is_void(Zi)){
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
  P = circularPupFunction(x,y,tel.obs);
  OTF_stats = autocorrelation(P*exp(1i*phi));

  //normalization
  surface_pup_m2 = tel.diam^2*(1-tel.obs^2)*pi/4;
  surface_pup_pix = surface_pup_m2 / tel.pixSize^2;
  factnorm = surface_pup_pix^2;
  OTF_stats /= factnorm;

  PSF_stats = fft(OTF_stats).re;
  
  return OTF_stats;
}

func computeFakeOTFncpa(SR_ncpa)
/* DOCUMENT
 */
{
  //defining a fake PSD
  N = tel.nPixels
  PSD_ncpa = array(0.0, N, N);

  WFE_ncpa = sqrt(- (cam.lambda*1e9/2/pi)^2 * log(SR_ncpa));
  fact = 2*pi* WFE_ncpa / tel.fourierPixSize / (cam.lambda*1e9);
  fact = fact^2;
  
  //defining the frequency mask of the DM
  k = computeSpatialFreqRad(N, tel.fourierPixSize, kx, ky);
  msk = defineDmFrequencyArea(k, kx, ky,"circle", tel.pitch );
  ndm = where( msk );
  northo = where( !msk );
  //filling the PSD
  tot = numberof(northo);      
  PSD_ncpa(northo) = fact / double(tot);

  //determining the phase structure function
  PSD_ncpa(N/2,N/2) = 0.0;
  PSD_ncpa(N/2,N/2) = -sum(PSD_ncpa);
  DPHI_ncpa = 2*abs(fft(roll(PSD_ncpa))) *  tel.fourierPixSize^2;

  //computing OTF
  OTF_tel = OTF_telescope(tel.diam,tel.obs,tel.nPixels,tel.pixSize);
  msk = OTF_tel/max(OTF_tel) > 1e-5;
  OTF_ncpa = exp(-0.5*DPHI_ncpa);

  //normalization
  OTF_ncpa *= SR_ncpa/sum(OTF_ncpa*OTF_tel)
    

  return OTF_ncpa;
}
func getOTFncpa(nPixels,procDir,&SR_bench,&PSF_ncpa,disp=)
{

  // .... finding the directories where ncpa calibration are storaged
  tmp = giveNCPADataDir(procDir);
  DirNCPA = tmp(1);
  timencpa = tmp(2);
  dataDirRoot = dataDir +  DirNCPA;
  // .... loading, cleaning (background extraction, dead pixel processing) and cropping the best PSF got on bench
  PSF_ncpa = processCanaryPSF(timencpa,SR_bench,box=nPixels,disp=disp);
  OTF_ncpa = fft(roll(PSF_ncpa)).re;

  return OTF_ncpa;

}

func giveNCPADataDir(Dir){

  if(strpart(Dir,1:7)=="2013_09"){
    DirNCPA = "2013_09_15_onsky/";
    suffirbench = "19h57m41s";
  }
  else if(strpart(Dir,1:7)=="2013_07"){
    DirNCPA = "2013_07_22_onsky/";
    suffirbench = "17h41m59s";
  }

  else if( strpart(Dir,1:7)=="2012_10" ){
    DirNCPA = "2012_10_04_onsky/";
    suffirbench = "20h43m28s";
  }

  
  return [DirNCPA,suffirbench];
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

func circularPupFunction(x,y,obs){
  r = abs(x,y);
  return r<=1.0 & r>obs;
}

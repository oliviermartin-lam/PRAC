
func PRAC_errorbreakdown(void,verb=)
/*DOCUMENT 

 */
{
                                                      

  // NCPA: value calibrated on bench.
  budget.ncpa   = 115;
  psf.SR_ncpa   = budget.SRncpa = var2sr(budget.ncpa,cam.lambda);

  slopesdis  = *rtc.slopes_dis;
  slopesdis -= slopesdis(,avg);

  //................. RESIDUAL ...................//
  sigTS_full = sigTS = 0;
  if(rtc.obsMode == "SCAO"){
    Son        = (*rtc.slopes_res)(slrange(rtc.its),);
    sigTS_full = computesResidual(Son,full=1);
    sigTS      = sqrt(sum(sigTS_full^2));
  }

  
  //................. TOMOGRAPHY .................//
  sigTomoIR_full = sigTomoIR = 0;
  if(rtc.obsMode == "MOAO"){
    covPara        = *covMatrix.parallel;
    sigTomoIR_full = computesTomoError(covPara, *rtc.R, rtc.its,full=1);
    sigTomoIR      = sqrt(sum(sigTomoIR_full^2));
    budget.tomo    = sigTomoIR;
  }
  //................. VED .....................//
if(rtc.obsMode == "MOAO"){
  for(l=1;l<=atm.nLayers;l++){
    covLearn_hl   = covMat1layer(rtc.nWfs,atm.cnh(l),atm.altitude(l),atm.l0h(l),sys.xshift,sys.yshift,sys.magnification,sys.theta,loworder=1);
    covLearn_hl   = handle_tilt_from_wfstype(covLearn_hl);
    budget.ved(l) = computesTomoError(covLearn_hl, *rtc.R, rtc.its);
  }
 }
  
  //................. ALIASING .................//
  covAlias = *covMatrix.aliasing;
  if(rtc.obsMode == "MOAO"){
    covAlias(slrange(rtc.its),) = 0;
    covAlias(,slrange(rtc.its)) = 0;
  }
  sigAliasIR_full = computesTomoError(covAlias, *rtc.R, rtc.its,full=1);
  sigAliasIR      = sqrt(sum(sigAliasIR_full^2));
  budget.alias    = sigAliasIR;
  
  //................. NOISE .................//
  //computes the time-filtering attenuation from the loop
  fgt = sqrt(rtc.loopGain/(2-rtc.loopGain));
 
  if(rtc.obsMode == "MOAO"){
    sigNoise        = computesNoisesError(takesDiag(*covMatrix.noise),*rtc.R,full=1);
    sigNoiseIR_full = fgt*sigNoise(,2);
  }else if(rtc.obsMode == "SCAO"){
    sigNoise        = computesNoisesError(takesDiag(*covMatrix.noiseCL),[],full=1);
    sigNoiseIR_full = fgt*sigNoise(,1);
  }
  sigNoiseIR = sqrt(sum(sigNoiseIR_full^2));
  budget.noise = sigNoiseIR;

  //................. BANDWIDTH .................//
  if(rtc.obsMode == "MOAO"){
    Son = (*rtc.R)(,+)*(slopesdis(norange(rtc.its),))(+,);  
    // Computing the servo lag error noise-free
    sigBW_full = computesBandwidthError(Son,rtc.obsMode,full=1);
    //subtracting aliasing contribution in OL
    sigBW_full = SQRT(sigBW_full^2 - sigAliasIR_full^2);
  }else{
    sigBW_full = SQRT(sigTS_full^2 - sigNoise^2 - sigAliasIR_full^2);
  }

   sigBW      = sqrt(sum(sigBW_full^2));
   budget.bw  = sigBW;
  
  //................. FITTING .................//
  
  sigFit = computesFittingError(tel.diam,atm.r0,atm.lambda*1e9);
  budget.fit = sigFit;

  //................. GO-TO ERROR .................//
  sigOL = 0;
  if(rtc.obsMode=="MOAO"){
    sigOL = computesOLError(*rtc.slopes_dis,*rtc.R);
  }
  budget.ol = sigOL;
  
  //................. STATIC .................//
  nm             = dimsof(*sys.slopesToZernikeMatrix)(2);
  sigStatic_full = array(0.,nm);
  sigStatic      = 0;
  if(rtc.obsMode == "MOAO"){
    sigStatic_full = computesStatic(*rtc.slopes_res,full=1);//in nm rms
    sigStatic      = sqrt(sum(sigStatic_full^2));
  }
  budget.static = sigStatic;
  
  //................. NCPA .................//                                           
  sigNCPA = budget.ncpa;

  //................. TOTAL .................//
  ir_rms2  = sigTomoIR^2 + sigAliasIR^2 +  sigNoiseIR^2 + sigBW^2 + sigFit^2 + sigOL^2;
  ir_total = ir_rms2 + sigStatic^2 + sigNCPA^2;

  budget.res = sqrt(ir_total);
  
  //adds the variances of the modes
  nzer   = numberof(sigBW_full);
  sigPHI =  sigBW_full^2 +  sigTomoIR_full^2 + sigAliasIR_full^2 +  sigNoiseIR_full^2 + sigStatic_full^2 + (sigOL^2)/nzer;
  sigPERP = sigFit^2 + sigNCPA^2;

  //computes the size of the halo
  r0tot    = atm.r0 *(cam.lambda/atm.lambda)^(1.2);
  sizeHalo = (1. + (tel.diam/r0tot)^2.);
  lamb2    = (2*pi/cam.lambda/1e9)^2;

  //................. SR .................//
  
  //Marechal approximation
  budget.SRmar = exp(-ir_total*lamb2);

  //Parenti approximation
  SRttr = exp(-(sum(sigPHI(3:)) + sigPERP)*lamb2);
  sigTT = (sum(sigPHI(1:2)))*lamb2;
  //expression of Parenti 1994
  budget.SRpar = SRttr/(1. + sigTT) + (1. - SRttr)/sizeHalo;

  //Born approximation
  SRperp        = exp(-sigPERP*lamb2);
  sigPara       = sum(sigPHI)*lamb2;
  budget.SRborn = SRperp/(1. + sigPara) + (1. - SRperp)/sizeHalo;

}




/*
 _   _  ___ ___ ____  _____ 
| \ | |/ _ \_ _/ ___|| ____|
|  \| | | | | |\___ \|  _|  
| |\  | |_| | | ___) | |___ 
|_| \_|\___/___|____/|_____|
                            
*/

func computesNoisesError(noiseVar,R,nho=,full=)
/* DOCUMENT [sigma_noisets,sigma_noiseR] = computeNoisesError(diag_noise,R)

  Computes noise error from the TS and the off-axis noise propagated
  through the reconstructor R in nm rms.
  
*/
{
  if(is_void(nho))
    nho = 3;

  Cnn            = array(0.,rtc.nSlopes,rtc.nSlopes);
  takesDiag(Cnn) = noiseVar;
    
  if(is_void(R)){
    mrz = *sys.slopesToZernikeMatrix;
    Cnn = Cnn(slrange(rtc.its),slrange(rtc.its));
    Czz = mrz(,+)*(Cnn(,+)*mrz(,+))(+,);
    nzz = dimsof(Czz)(0);
    // units ...
    units = 0.5*tel.diam*1e9/radian2arcsec;
    diagonale = takesDiag(Czz)*units^2;
    if(full){
      return SQRT(diagonale);
    }else{
      return SQRT((sum(diagonale)));
    }

  }else{
    //computes TS noise error
    tsnoise        = computesTomoError(Cnn,0*R,full=full);
    //Computes Tomo noise error
    Cnn(slrange(rtc.its),) = 0;
    Cnn(,slrange(rtc.its)) = 0;
    tomonoise      = computesTomoError(Cnn,R,full=full); // in nm rms
    if(!full){
      sigNoiseTS   = tsnoise;
      sigNoiseTomo = tomonoise(1);
      return [sigNoiseTS,sigNoiseTomo];
    }else
      return [tsnoise,tomonoise];
  }
}



 /*
 ____    _    _   _ ______        _____ ____ _____ _   _ 
| __ )  / \  | \ | |  _ \ \      / /_ _|  _ \_   _| | | |
|  _ \ / _ \ |  \| | | | \ \ /\ / / | || | | || | | |_| |
| |_) / ___ \| |\  | |_| |\ V  V /  | || |_| || | |  _  |
|____/_/   \_\_| \_|____/  \_/\_/  |___|____/ |_| |_| |_|
*/

func computesBandwidthError(Son, obsmode,nho=,full=)
/*DOCUMENT BWerr = computesBandwidthError(dataset, retard, gain, icam,obsmode,verbose=)

  Estimates the Bandwidth error on the icam part of dataset from a
  simulation of the MOAO loop. The BW error is deduced by the
  residual of the simulation ( BW error and noise error) minus the
  noise error ( measured and propagated through the loop )

  SEE ALSO: filteredTomoNoiseSeenByTS_moao
 */
{

  if(is_void(nho))
    nho = 3;

 

  //Simulating the system
  if(obsmode == "MOAO"){
    sp    = (1-rtc.loopGain)*roll(Son,[0,1]) + rtc.loopGain*Son;
  
  
    a     = rtc.delay*rtc.Fe ;
    r     = Son - (a*roll(sp,[0,2]) + (1-a)*roll(sp,[0,1]));
    r     = r(,4:);

    //Expansion over Zernike polynomials
    mrz   = *sys.slopesToZernikeMatrix;
    units = 0.5*(tel.diam)*1e9/radian2arcsec;
    rz    = (mrz(,+)*r(+,))*units;
  

    //Computing noise Zernike expansion in nm^2
    nz = getNoiseZernike(Son,rtc.its,full=1,arc=1);
    fact  = filteringNoiseFactor(rtc.loopGain,rtc.delay*rtc.Fe,rtc.obsMode);
  
    //subtracting time-filtered noise variance
    rz    = rz(,rms)^2 - fact*nz;

    if(!full){
      //estimation of the residual in nm
      sigBW = sqrt(rz(sum));
      sigBW_TTR = sqrt(sum(rz(nho:)));
      return [sigBW,sigBW_TTR];
    }else{
      return sqrt(rz);
    }
  }


  /*
  else{

  
  //Temporal frequencies domain
  nFrames  = rtc.nFrames;
  nu       = (span(0,1,nFrames)*rtc.Fe);
  dnu      = avg(nu(dif));
 
  //...... Computes temporal DSP for each Zernike modes .....//
  Son     -= Son(,avg);
  mrz      = *sys.slopesToZernikeMatrix;
  //derives noise and Zernike basis
  units    = 0.5*(tel.diam)*1e9/radian2arcsec;
  data_z   = mrz(,+) * Son(+,) * units;

  //Temporal dsp of the Zernike modes of the noise
  zDSP     = abs(fft(data_z ,[0,1]))^2. /nFrames^2;

  //transfer function
  if(obsmode == "MOAO"){
    hcor = hcorMoao(nu,rtc.Fe,rtc.delay,rtc.loopGain,rtc.BP);
  }else if(obsmode =="SCAO"){
    hcor = hcorScao(nu,rtc.Fe,rtc.delay,rtc.loopGain,rtc.BP);
  }

  //temporal dsp of the Zernike modes of the residual turbulence 
  zDSP = (abs(hcor(-,))^2 * zDSP);
  zDSP = zDSP(,1:nFrames/2);
  zDSP = simpsonIntegral(zDSP,dnu);


  if(!full){
    //estimation of the residual in nm
    sigBW = sqrt(zDSP(sum));
    sigBW_TTR = sqrt(sum(zDSP(nho:,)));
  return [sigBW,sigBW_TTR];
  }else
    return sqrt(zDSP);
  }
  */
  }



/*
 _____ ___ _____ _____ ___ _   _  ____ 
|  ___|_ _|_   _|_   _|_ _| \ | |/ ___|
| |_   | |  | |   | |  | ||  \| | |  _ 
|  _|  | |  | |   | |  | || |\  | |_| |
|_|   |___| |_|   |_| |___|_| \_|\____|
*/

func computesFittingError(D,rzero,lambda)
/*DOCUMENT sigma_fitting = computesFittingError(Dpup,rzero,lambdaWFS)

  Returns the fitting error in nm rms at the wavelength lambdaWFS in nm with Dpup
  the diameter of the pupil in meters and rzero the Fried parameters
  in meters
 */
{
  // sqrt(0.257*(35.^(-5./6.))/2/pi = 0.01834108434934971638
  return sqrt(0.0122) * (D/rzero)^(5/6.) * lambda/2./pi;
  //return sqrt(0.257*(35.^(-5./6.)))/2/pi*lambda*(D / rzero)^(5/6.);
}

 /* 
  ____  ___      _____ ___        _____ ____  ____   ___  ____  
 / ___|/ _ \    |_   _/ _ \      | ____|  _ \|  _ \ / _ \|  _ \ 
| |  _| | | |_____| || | | |_____|  _| | |_) | |_) | | | | |_) |
| |_| | |_| |_____| || |_| |_____| |___|  _ <|  _ <| |_| |  _ < 
 \____|\___/      |_| \___/      |_____|_| \_\_| \_\\___/|_| \_\
                                                                
  */
func computesOLError(dataset,R)
/*DOCUMENT
 */
{
  mrz = *sys.slopesToZernikeMatrix;
  units = 0.5*(tel.diam)*1e9/radian2arcsec;
      
  //determines the tomographic prediction of the TS measurements
  TSest = R(,+)*(dataset(norange(rtc.its),))(+,);
  //Derives on Zernike
  TSest_z = mrz(,+)*TSest(+,) *units;
  //takes variance
  sigma_WF = sqrt(sum(TSest_z(,rms)^2));
  //Determines OL error in taking 5% of the wavefront error
  sigOL = 0.05*sigma_WF;

  return sigOL;

}

 /*
  ____ _____  _  _____ ___ ____ 
/ ___|_   _|/ \|_   _|_ _/ ___|
\___ \ | | / _ \ | |  | | |    
 ___) || |/ ___ \| |  | | |___ 
|____/ |_/_/   \_\_| |___\____|
                               
  */
func computesStatic(dataset,nho=,full=)
/* DOCUMENT [TSrms,TSrms_wTT,TSstatic,TSstatic_wTT] = computesTSresiduals(dataset)

   Determines the rms and static residual od the TS measurements in
   dataset.
   
 */
{
  if(is_void(nho))
    nho=3;
  icam = rtc.its;
  //takes mrz
  mrz = *sys.slopesToZernikeMatrix;
  units = 0.5*(tel.diam)*1e9/radian2arcsec;
  ZZ = mrz(,+)*(dataset(slrange(icam),))(+,) *units;

  //Takes dynamic and static residuals
  sigStatic_full = ZZ(,avg);   // static abs, in nm rms

  if(!full){
    //Separate TT and high order components
    sigStatic = (ztstat^2)(sum)^0.5;    // Total static error
    sigStatic_TTR = (ztstat(3:)^2)(sum)^0.5;
    sigStatic_HO = (ztstat(nho:)^2)(sum)^0.5;    // Static error without Tip Tilt

    return [sigStatic, sigStatic_TTR, sigStatic_HO_wTT];
  }else
    return sigStatic_full;

}

/*
 ____  _____ ____ ___ ____  _   _   _    _     
|  _ \| ____/ ___|_ _|  _ \| | | | / \  | |    
| |_) |  _| \___ \| || | | | | | |/ _ \ | |    
|  _ <| |___ ___) | || |_| | |_| / ___ \| |___ 
|_| \_\_____|____/___|____/ \___/_/   \_\_____|
                                               
*/

func computesResidual(Son,nho=,full=)
/* DOCUMENT [TSrms,TSrms_wTT] = computeResidual(dataset)

   Expandes dataset over Zernike polynomials to get
   the rms value in nm.
   
 */
{
  if(is_void(nho))
    nho=3;
  
  // Mean removal
  Son -= Son(,avg);
  
  //Expansion over Zernike polynomials
  mrz   = *sys.slopesToZernikeMatrix;
  units = 0.5*(tel.diam)*1e9/radian2arcsec;
  ZZ    = mrz(,+)*Son(+,) *units;

  //Getting the rms value
  sigRes_full = ZZ(,rms)^2;  

  if(!full){
    //Separate TT and high order components
    sigRes     = sqrt(sigRes_full(sum));
    sigRes_TTR = sqrt((sigRes_full(3:))(sum));
    sigRes_HO  = sqrt((sigRes_full(1:2))(sum))

    return [sigStatic, sigStatic_TTR, sigStatic_HO_wTT];
  }else
    return sqrt(sigRes_full);

}

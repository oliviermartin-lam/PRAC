
func PRAC_errorbreakdown(void,verb=)
/*DOCUMENT 

 */
{
                                                      

  slopesdis = *data.slopes_dis;
  slopesdis -= slopesdis(,avg);
  
  //................. TOMOGRAPHY .................//
  sigTomoIR_full = sigTomoIR = 0;
  if(data.rtc.obsmode == "MOAO"){
    covPara = *data.covmatrix.covPara - *data.covmatrix.covTracking;
    sigTomoIR_full = computesTomoError(covPara, *data.rtc.R, data.its,full=1);
    sigTomoIR = sqrt(sum(sigTomoIR_full^2));
    data.budget.tomo = sigTomoIR;
  }
  //................. VED .....................//
if(data.rtc.obsmode == "MOAO"){
  for(l=1;l<=data.learn.nl;l++){
    covLearn_hl = covMat1layer(data.nwfs,l,data.learn.cnh(l),data.learn.altitude(l),data.learn.l0h(l),loworder=1);
    covLearn_hl = handle_tilt_from_wfstype(covLearn_hl);
    data.budget.ved(l) = computesTomoError(covLearn_hl, *data.rtc.R, data.its);
  }
 }
  
  //................. ALIASING .................//
  covAlias = *data.covmatrix.covAlias;
  if(data.rtc.obsmode == "MOAO"){
    covAlias(slrange(data.its),) = 0;
    covAlias(,slrange(data.its)) = 0;
  }
  sigAliasIR_full = computesTomoError(covAlias, *data.rtc.R, data.its,full=1);
  sigAliasIR = sqrt(sum(sigAliasIR_full^2));
  data.budget.alias = sigAliasIR;
  
  //................. NOISE .................//
  //computes the time-filtering attenuation from the loop
  fgt = sqrt(filteringNoiseFactor(data.rtc.gain,data.rtc.delay*data.rtc.Fe,data.rtc.obsmode));
  if(data.rtc.obsmode=="MOAO"){
    sigNoiseIR_full =  fgt*(computesNoisesError(*data.turbu.varNoise,*data.rtc.R,full=1))(,2);
  }else if(data.rtc.obsmode =="SCAO"){
    sigNoiseIR_full = getNoiseZernike(slopesdis(slrange(data.its),),data.its,full=1,arc=1)*1e6;
    sigNoiseIR_full = fgt*sqrt(sigNoiseIR_full);
  }
  sigNoiseIR = sqrt(sum(sigNoiseIR_full^2));
  data.budget.noise = sigNoiseIR;

  //................. BANDWIDTH .................//
  if(data.rtc.obsmode == "MOAO"){
    Son = (*data.rtc.R)(,+)*(slopesdis(norange(data.its),))(+,);
  }else if (data.rtc.obsmode =="SCAO"){
    Son = slopesdis(slrange(data.its),);
  }
  
  sigBW_full = SQRT(computesBandwidthError(Son,data.rtc.obsmode,full=1)^2 - sigNoiseIR_full^2);
  sigBW = sqrt(sum(sigBW_full^2));
  data.budget.bw = sigBW;

  
  //................. FITTING .................//
  sigFit = computesFittingError(data.tel.diam,data.turbu.r0vis,data.lambda_vis*1e9);
  data.budget.fit = sigFit;

  //................. GO-TO ERROR .................//
  sigOL = 0;
  if(data.rtc.obsmode=="MOAO"){
    sigOL = computesOLError(*data.slopes_dis,*data.rtc.R);
  }
  data.budget.ol = sigOL;
  
  //................. STATIC .................//
  nm = dimsof(*data.wfs(data.its).mrz)(2);
  sigStatic_full = array(0.,nm);
  sigStatic = 0;
  if(data.rtc.obsmode == "MOAO"){
    sigStatic_full = computesStatic(*data.slopes_res,full=1);//in nm rms
    sigStatic = sqrt(sum(sigStatic_full(3:)^2));
  }
  data.budget.static = sigStatic;
  
  //................. NCPA .................//                                           
  sigNCPA = data.budget.ncpa;

  //................. TOTAL .................//
  ir_rms2 =  sigTomoIR^2 + sigAliasIR^2 +  sigNoiseIR^2 + sigBW^2 + sigFit^2 + sigOL^2;
  ir_total = ir_rms2 + sigStatic^2 + sigNCPA^2;

  data.budget.res = sqrt(ir_total);
  
  //adds the variances of the modes
  sigPHI =  sigBW_full^2 +  sigTomoIR_full^2 + sigAliasIR_full^2 +  sigNoiseIR_full^2 ;
  //adds OL and static
  sigPHI(3:) += sigStatic_full(3:)^2 + (sigOL^2)/numberof(sigPHI(3:));
  sigPERP = sigFit^2 + sigNCPA^2;

  //computes the size of the halo
  sizeHalo = (1. + (data.tel.diam/data.turbu.r0ir)^2.);
  lamb2 = (2*pi/data.camir.lambda_ir/1e9)^2;

  //................. SR .................//
  
  //Marechal approximation
  data.budget.SRmar = 100*(exp(-ir_total*lamb2));

  //Parenti approximation
  SRttr = exp(-(sum(sigPHI(3:)) + sigPERP)*lamb2);
  sigTT = (sum(sigPHI(1:2)))*lamb2;
  //expression of Parenti 1994
  data.budget.SRpar = 100*(SRttr/(1. + sigTT) + (1. - SRttr)/sizeHalo);

  //Born approximation
  SRperp = exp(-sigPERP*lamb2);
  sigPara = sum(sigPHI)*lamb2;
  data.budget.SRborn = 100*(SRperp/(1. + sigPara) + (1. - SRperp)/sizeHalo);

}









/*
 _   _  ___ ___ ____  _____ 
| \ | |/ _ \_ _/ ___|| ____|
|  \| | | | | |\___ \|  _|  
| |\  | |_| | | ___) | |___ 
|_| \_|\___/___|____/|_____|
                            
*/

func computesNoisesError(noiseVar,R,nho=,full=)
/*DOCUMENT [sigma_noisets,sigma_noiseR] = computeNoisesError(diag_noise,R)

  Computes noise error from the TS and the off-axis noise propagated
  through the reconstructor R.
  
*/
{
  if(is_void(nho))
    nho = 3;
  //defines noises covariance matrices
  Cnn = array(0.,data.Nslopes,data.Nslopes);
  takesDiag(Cnn) = noiseVar;
  //computes TS noise error
  tsnoise = computesTomoError(Cnn,0*R,full=full);
  //Computes Tomo noise error
  Cnn(slrange(data.its),) = 0;
  Cnn(,slrange(data.its)) = 0;
  tomonoise = computesTomoError(Cnn,R,full=full); // in nm rms
  if(!full){
    sigNoiseTS = tsnoise(1);
    sigNoiseTS_TTR = tsnoise(2);
    
    sigNoiseTomo = tomonoise(1);
    sigNoiseTomo_TTR = tomonoise(2);
    
    return [sigNoiseTS,sigNoiseTS_TTR,sigNoiseTomo,sigNoiseTomo_TTR];
  }else
    return [tsnoise,tomonoise];
}



 /*
 ____    _    _   _ ______        _____ ____ _____ _   _ 
| __ )  / \  | \ | |  _ \ \      / /_ _|  _ \_   _| | | |
|  _ \ / _ \ |  \| | | | \ \ /\ / / | || | | || | | |_| |
| |_) / ___ \| |\  | |_| |\ V  V /  | || |_| || | |  _  |
|____/_/   \_\_| \_|____/  \_/\_/  |___|____/ |_| |_| |_|
*/

func computesBandwidthError(dataset, obsmode,nho=,full=)
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
  
  //Temporal frequencies domain
  Nframes = dimsof(dataset)(0);
  nu = data.rtc.Fe*(indgen(Nframes)-1)/Nframes;
  dataset -= dataset(,avg);
  //...... Computes temporal DSP for each Zernike modes .....//
  mrz = *data.wfs(data.its).mrz;
  //derives noise and Zernike basis
  units = ((1/206265.)*(0.5*(data.tel.diam)*1e9));
  data_z = mrz(,+) * dataset(+,) * units;
    
  //Temporal dsp of the Zernike modes of the noise
  dspdata_z = fft(data_z ,[0,1]);

  //transfer function
  if(obsmode == "MOAO"){
    hcor = hcorMoao(nu,data.rtc.Fe,data.rtc.delay,data.rtc.gain,data.rtc.BP);
    dspdatares_z = 4*(abs(hcor(-,)*dspdata_z)^2./Nframes^2.)(,1:Nframes/2);
  }else if(obsmode =="SCAO"){
    hcor = hcorScao(nu,data.rtc.Fe,data.rtc.delay,data.rtc.gain,data.rtc.BP);
    dspdatares_z = (abs(hcor(-,)*dspdata_z)^2./Nframes^2.)(,1:Nframes/2);
  }
  //temporal dsp of the Zernike modes of the residual turbulence
  

  if(!full){
    //estimation of the residual in nm
    sigBW = sqrt(sum(sigCross_full) - sum(dspdatares_z));
    sigBW_TTR = sqrt(sum(dspdatares_z(nho:,)));
  return [sigBW,sigBW_TTR];
  }else
    return sqrt(dspdatares_z(,sum));
}



/*
 _____ ___ _____ _____ ___ _   _  ____ 
|  ___|_ _|_   _|_   _|_ _| \ | |/ ___|
| |_   | |  | |   | |  | ||  \| | |  _ 
|  _|  | |  | |   | |  | || |\  | |_| |
|_|   |___| |_|   |_| |___|_| \_|\____|
*/  
func computesFittingError(Dpup,rzero,lambdaWFS)
/*DOCUMENT sigma_fitting = computesFittingError(Dpup,rzero,lambdaWFS)

  Returns the fitting error in nm rms at the wavelength lambdaWFS in nm with Dpup
  the diameter of the pupil in meters and rzero the Fried parameters
  in meters
 */
{
  // sqrt(0.257*(35.^(-5./6.))/2/pi = 0.01834108434934971638
  return 0.01834108434934971638*lambdaWFS*(Dpup / rzero)^(5/6.);
}

 /* 
  ____  ___      _____ ___        _____ ____  ____   ___  ____  
 / ___|/ _ \    |_   _/ _ \      | ____|  _ \|  _ \ / _ \|  _ \ 
| |  _| | | |_____| || | | |_____|  _| | |_) | |_) | | | | |_) |
| |_| | |_| |_____| || |_| |_____| |___|  _ <|  _ <| |_| |  _ < 
 \____|\___/      |_| \___/      |_____|_| \_\_| \_\\___/|_| \_\
                                                                
  */
func computesOLError(dataset,R,pix=,lib=)
/*DOCUMENT
 */
{
  mrz = *data.wfs(data.its).mrz;
  units = ((1/206265.)*(0.5*(data.tel.diam)*1e9));
      
  //determines the tomographic prediction of the TS measurements
  TSest = R(,+)*(dataset(norange(data.its),))(+,);
  //Derives on Zernike
  TSest_z = mrz(,+)*TSest(+,) *units;
  //takes variance of the tip-tilt wavefront only
  sigma_TT = sqrt(sum(TSest_z(3:,rms)^2));
  //Determines OL error in taking 5% of the wavefront error on high order only
  sigOL = 0.05*sigma_TT;

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
  icam = data.its;
  //takes mrz
  mrz = *data.wfs(icam).mrz;
  //Derives WFS measurements to Zernike basis
  units = ((1/206265.)*(0.5*(data.tel.diam)*1e9));
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

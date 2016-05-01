pathPRAC = "/home/omartin/CANARY/Algorithms/PRAC/";

include, pathPRAC + "pracDm.i"; 
include, pathPRAC + "pracPSD.i";
include, pathPRAC + "fourierAdaptiveOptics.i";
include, pathPRAC + "pracTomography.i";
include, pathPRAC + "pracLmfit.i";
include, pathPRAC + "pracLearn.i";
include, pathPRAC + "pracErrorbreakdown.i";
include, pathPRAC + "pracAtmosphere.i";
include, pathPRAC + "pracDirectPsfr.i";
include, pathPRAC + "pracMakeOTF.i";

include, pathPRAC + "canaryUtils.i";
include, pathPRAC + "telescopeUtils.i"; 
include, pathPRAC + "noiseUtils.i";
include, pathPRAC + "imageUtils.i";
include, pathPRAC + "contextUtils.i";
include, pathPRAC + "fitsUtils.i";
include, pathPRAC + "mathUtils.i";
include, pathPRAC + "displayUtils.i";
include, pathPRAC + "constantUtils.i";
include, pathPRAC + "zernikeUtils.i";



func pracMain(timedata,Dir=,psfr=,averageMode=,budgetonly=,verb=,disp=,writeRes=)
/* DOCUMENT res = pracMain(timedata,Dir=,psfr=,averageMode=,verb=,disp=,writeRes=)

  Launches the PSF reconstruction for processing the CANARY on-sky data slopestl acquired at "timedata"
  and storaged into the directory specified by the keyword "Dir"
  using the method given by the keyword "psfr."

  INPUTS:

  - timedata (string)   : local hour of the slopestl file
  - Dir (string)        : Night directory in which the slopestl are storaged.
  It depends on how you've nammed your own data directories.
  - psfr (string)       : Reconstruction method: "estimation" for analytic reconstruction,
  "rtc" for the RTC-based method and "ts" for the TS-based method.
  Those three methods are available in SCAO and MOAO (any reconstructors).
  - verb (boolean)      : set to 1 to get information about what the algorithm is doing and final results
  - disp (boolean)      : set to 1 to get graphic results
  - writeRes (boolean)  : set to 1 to save automatically the returning pointer res;

  OUTPUTS:

  - res                 : this is the returning pointer of 12 elements:
     
  res(1) = &strchar([timedata,rtc.aoMode,rtc.obsMode,rtc.recType]);
  res(2) = &[atm.r0,atm.L0,atm.v,tf];
  res(3) = &[atm.cnh,atm.altitude,atm.l0h,atm.vh];
  res(4) = &[sys.tracking];
  res(5) = &[wfs.x,wfs.y,sys.xshift,sys.yshift,sys.magnification,sys.theta,sys.centroidGain];
  res(6) = &[*psf.sky,*psf.res,*psf.diff,*psf.ncpa];
  res(7) = &[*psf.EE_sky,*psf.EE_res];
  res(8) = &[*otf.sky,*otf.res];
  res(9) = &[budget.res,budget.fit,budget.bw,budget.tomo,budget.noise,
             budget.alias,budget.static,budget.ncpa,budget.ol];
  res(10) = &budget.ved;
  res(11) = &([psf.SR_sky,psf.SR_res,budget.SRmar,budget.SRpar,budget.SRborn,
                  psf.SR_tomo,psf.SR_fit,psf.SR_bw,psf.SR_stats,psf.SR_ncpa,
                  psf.FWHM_sky,psf.FWHM_res,psf.chi2]);
  res(12) = &[psf.skyMoffatProf,psf.resMoffatProf];

  EXAMPLES:
  
  res = pracMain("00h15m36s",Dir="2013_09_13_onsky/",verb=1,disp=1,psfr = "estimation")

  performs the PSF-R using the estimation method on the slopestl file acquired at 00h15m36s 
  and storaged into the night directory "2013_09_13_onsky/".

  SEE ALSO: concatenatePracResults()
 */
{
  
  if(psfr!="estimation" & psfr!="ts" &psfr!="rtc" & psfr!="all" & !budgetonly){
    write,"Please make your choice between estimation, TS and RTC PSF-R methods";
    return [];
  }
  
  // Directory management
  include, "pracConfig.i",1;
  if(Dir){
    procDir = Dir;
    dataDirRoot = dataDir + procDir;
  }
  if(is_void(averageMode))
    averageMode = "Vii";
  
  tic,1;
  
  /////////////////////
  //.... Initializing the data struct
  //////////////////////////////////////////////////////////
  tic,2;
  include, "pracStructConfig.i",1;
  define_structs,timedata,verb=verb;
  t_init = tac(2);
  
  /////////////////////
  // .... Error breakdown computation
  //////////////////////////////////////////////////////////
  tic,3;
  PRAC_errorbreakdown,verb=verb;
  t_budg = tac(3);
  if(budgetonly)
    return concatenatePracResults();

  /////////////////////
  // .... Perfect telescope OTF
  //////////////////////////////////////////////////////////
  tic,4;
  telOTF = pracTelescopeOTF(tel.obs,tel.nPixels);
  otf.tel = &roll(telOTF);

  /////////////////////
  // .... OTF from static aberrations
  //////////////////////////////////////////////////////////

  //Common path static aberrations without telescope
  statOTF      = staticOTF(statPSF,statSR);
  otf.static   = &roll(statOTF);
  psf.SR_stats = statSR;

  /////////////////////
  // .... OTF from fitting error
  //////////////////////////////////////////////////////////
   
  fitOTF      = roll(PSD2OTF(fittingPSD(fao),tel.fourierPixSize));
  otf.fit     = &roll(fitOTF);
  psf.SR_fit  = sum((*otf.fit) * (*otf.tel))/sum(*otf.tel);

  /////////////////////
  // .... OTF from bench-calibrated NCPA
  //////////////////////////////////////////////////////////
    
  //Getting the best bench PSF
  ncam     = cam.nPixelsCropped;
  ncpaOTF  = getOTFncpa(ncam,procDir,SRbench,PSF_ncpa,disp=disp);
  ncpaOTF  = roll(ncpaOTF);
  //ncpaOTF  = roll(*otf.tel);
  otf.ncpa = &roll(ncpaOTF);
  psf.ncpa = &PSF_ncpa; 

  
  /////////////////////
  // .... PSF Reconstruction
  //////////////////////////////////////////////////////////
  
  if(psfr == "estimation"){

    resOTF = runEstimationMethod(psf,otf,rtc,averageMode=averageMode);

  }else if(psfr == "ts"){

    resOTF = runTSbasedMethod(psf,otf,rtc,covMatrix,averageMode=averageMode);

  }else if(psfr == "rtc"){

    resOTF = runRTCbasedMethod(psf,otf,rtc,covMatrix,averageMode=averageMode);

  }else if(psfr == "all"){
    
    resOTF_est = runEstimationMethod(psf,otf,rtc,averageMode=averageMode);
    resOTF_ts  = runTSbasedMethod(psf,otf,rtc,covMatrix,averageMode=averageMode);
    resOTF_rtc = runRTCbasedMethod(psf,otf,rtc,covMatrix,averageMode=averageMode);
  }
  otf.ts = &roll(resOTF);
  t_rec = tac(4);
  
  /////////////////////
  // .... Processing the sky PSF
  //////////////////////////////////////////////////////////
  tic,5;
  //to set the maximum intensity at the middle pixel
  skyOTF  = roll(*otf.sky);
  // Normalazing the PSF to the Strehl ratio
  skyPSF  = roll(fft(skyOTF).re);
  skyPSF /= sum(skyPSF);
  skyPSF /= tel.airyPeak;
  //Storaging information
  psf.sky      = &skyPSF;
  psf.SR_sky   = max(skyPSF);
  budget.SRsky = max(skyPSF);
  //Getting the FWHM
  psf.FWHM_sky    = getPsfFwhm(*psf.sky,cam.pixSize,asky,fit=2);
  psf.moffat_sky = asky;
  
  /////////////////////
  // ..... Reconstructed PSF processing
  //////////////////////////////////////////////////////////

  if(psfr == "all"){
    
    processReconstructedOTF,resOTF_ts ,telOTF,ncpaOTF,ats;
    SRts = psf.SR_res;FWHMts = psf.FWHM_res;

    processReconstructedOTF,resOTF_rtc,telOTF,ncpaOTF,artc;
    SRrtc = psf.SR_res;FWHMrtc = psf.FWHM_res;

    processReconstructedOTF,resOTF_est,telOTF,ncpaOTF,aest;
    SRest = psf.SR_res;FWHMest = psf.FWHM_res;

    pracRes = array(pointer,6);
    pracRes(1) = &[psf.SR_sky,SRrtc,SRest,SRts];
    pracRes(2) = &[psf.FWHM_sky,FWHMrtc,FWHMest,FWHMts];
    pracRes(3) = &asky;
    pracRes(4) = &artc;
    pracRes(5) = &aest;
    pracRes(6) = &ats;
    return pracRes;
  }
  
  processReconstructedOTF,resOTF,telOTF,ncpaOTF;
  /////////////////////
  // ..... Reconstruction residue
  //////////////////////////////////////////////////////////
    
  psf.diff    = &(*psf.sky - *psf.res);
  w           = where(*psf.sky >0);
  psf.diffAvg = ((*psf.diff)(w)/ (*psf.sky)(w))(avg);
  psf.diffRms = ((*psf.diff)(w) /(*psf.sky)(w))(rms);
  psf.chi2    = 100*sum(((*psf.diff)(w)))/sum((*psf.sky)(w));
  
   
  /////////////////////
  // ..... Getting Ensquared Energy and FWHM on both reconstructed/sky PSF
  //////////////////////////////////////////////////////////
  
  boxsize = EE_res = EE_sky = span(1,cam.nPixelsCropped-3.,cam.nPixelsCropped) * cam.pixSize;
  for(i=1; i<=cam.nPixelsCropped; i++){
    EE_res(i) = getEE(100*(*psf.res)/sum(*psf.res), cam.pixSize, boxsize(i));
    EE_sky(i) = getEE(100*(*psf.sky)/sum(*psf.sky), cam.pixSize, boxsize(i));
  }

  psf.EE_sky        = &EE_sky;
  psf.EE_res        = &EE_res;
  t_proc = tac(5);
 
  
  /////////////////////
  // .... Display and verbose
  //////////////////////////////////////////////////////////

  if(verb){
    tf = tac(1);
    printPracResults,psf,budget,t_init,t_budg,t_rec,t_proc,tf;
  }
  if(disp){
    displayPracResults,psf,otf,tel,atm,cam,budget;
  }

  /////////////////////
  // .... Concatenate results
  //////////////////////////////////////////////////////////
  
  pracResults = concatenatePracResults(psfr + "_" + averageMode,writeRes=writeRes);

  
  return pracResults;
}


func processReconstructedOTF(resOTF,telOTF,ncpaOTF,&ares)
/* DOCUMENT processReconstructedOTF,resOTF,telOTF,ncpaOTF

 */
{
  extern tel,cam,psf,otf;

  //Interpolating to the camera resolution
  resOTF  = roll(changeResolution(roll(resOTF),tel.pixSize,cam.pixSize,ncam*cam.pixSize));
 
  telOTF  = roll(changeResolution(roll(telOTF),tel.pixSize,cam.pixSize,ncam*cam.pixSize));
    
  
  //Handling the NCPA OTF acquired at 1550 nm to 1650 nm.
  ncpaPixSize = radian2arcsec * 1.55e-6 * 0.408472/4.2;
  ncpaOTF = roll(changeResolution(roll(ncpaOTF),ncpaPixSize,cam.pixSize,ncam*cam.pixSize));
  ncpaOTF /= sum(ncpaOTF)/(psf.SR_ncpa*sum(telOTF));

  //Multiplying by NCPA
  resOTF *= ncpaOTF;

  otf.res = &roll(resOTF);
  //grabbing the PSF
  PSF_res        = roll(fft(resOTF).re);  
  PSF_res        = rotate2(PSF_res,90);
  psf.SR_res     = sum(resOTF)/sum(telOTF);
  psf.res        = &(PSF_res/(max(PSF_res)/psf.SR_res));
  psf.FWHM_res   = getPsfFwhm(*psf.res,tel.pixSize,ares,fit=2);
  psf.moffat_res = ares;
  
}


func runEstimationMethod(psf,otf,rtc,averageMode=)
/* DOCUMENT runEstimationMethod(psf,otf,rtc)

 */
{

  /////////////////////
  // .... OTF from response delay time of the system
  //////////////////////////////////////////////////////////
  
  //servoOTF   = roll(servoLagOTF(rtc.obsMode));
  servoPSD   = servoLagPSD(fao);
  servoPSD   = enlargeSupport(servoPSD,tel.nTimes);
  servoOTF   = roll(PSD2OTF(servoPSD,tel.fourierPixSize));
  otf.bw     = &roll(servoOTF);
  psf.SR_bw  = sum((*otf.bw) * (*otf.tel))/sum(*otf.tel);
   
    
  /////////////////////
  // .... OTF from tomographic residue, mode = "Uij", "Vii" or "intersample"
  //////////////////////////////////////////////////////////////////////////////

  Cee = computesCeeMatrix(rtc.obsMode,para=para,verb=verb);
    
  //computes the covariance matrix of voltages in volts^2
  Cvv = propagateError2Dm(Cee,*rtc.mcsky );
        
  //Computing the OTF
  tomoOTF     = makeOTF(Cvv,dphi,averageMode=averageMode,verb=verb);
  otf.tomo     = &roll(tomoOTF);
  psf.SR_tomo  = sum((*otf.tomo) * (*otf.tel))/sum(*otf.tel);
    
  /////////////////////
  // .... OTF at the science camera location
  //////////////////////////////////////////////////////////
  
  resOTF =  statOTF * fitOTF * servoOTF * tomoOTF;
    
  return resOTF;
}


func runRTCbasedMethod(psf,otf,rtc,covMatrix,averageMode=)
/* DOCUMENT

 */
{

  if(rtc.obsMode == "MOAO"){
      
      /////////////////////////////////
      // .... Computing the best MMSE tomographic reconstructor
      /////////////////////////////////////////////////////

      //grabbing matrices
      Coffoff  = (*covMatrix.parallel)(norange(rtc.its),norange(rtc.its));
      Conoff   = (*covMatrix.parallel)(slrange(rtc.its),norange(rtc.its));
      Conon    = (*covMatrix.parallel)(slrange(rtc.its),slrange(rtc.its));
      R        = *covMatrix.R;
      
       /////////////////////////////////
      // .... MMSE reconstructor + Covariance of the error
      /////////////////////////////////////////////////////

      tmp = R(,+) * Conoff(,+);
      Cee = Conon - tmp - transpose(tmp) + R(,+)*(Coffoff(,+)*R(,+))(+,);
      
      /////////////////////////////////
      // .... Residual TS slopes from MMSE estimation
      /////////////////////////////////////////////////////
      
      // MMSE reconstruction
      soff = (*rtc.slopes_res)(norange(rtc.its),);
      son  = R(,+) * soff(+,);
      son  = mirror_SH7(son,1);
      mc   = *rtc.mcsky;
       
      // determining the delayed voltages
      V    = *rtc.volts;
      res  = mc(,+) * son(+,) + V;
      res -= res(,avg);
      Cres = res(,+) * res(,+)/dimsof(res)(0);
      Cres+= propagateError2Dm(Cee,mc);
      
      /////////////////////////////////
      // .... Determining the residual OTF
      /////////////////////////////////////////////////////
      
      resOTF = makeOTF(Cres,averageMode=averageMode,verb=verb);
           
    }else if(rtc.obsMode == "SCAO"){
           
      //Grabbing the voltages vectors in volts
      V    = (*rtc.volts)(,dif);
      // noise matrix
      Cnn  = (*covMatrix.noiseCL)(slrange(rtc.its),slrange(rtc.its));
      // aliasing matrix
      Caa  = (*covMatrix.aliasing)(slrange(rtc.its),slrange(rtc.its));
      //Computing the volts covariance matrix
      V   -= V(,avg);
      Cvv  = V(,+) * V(,+)/dimsof(V)(0);
      Cvv  = Cvv + propagateError2Dm(Cnn + Caa, *rtc.mcsky);     
      //Computing the OTF
      resOTF = makeOTF(Cvv,averageMode=averageMode,verb=verb);   
    }
    
    //Final OTF
    resOTF *= roll(*otf.fit) * roll(*otf.static);

    return resOTF;
}



func runTSbasedMethod(psf,otf,rtc,covMatrix,averageMode=)
/* DOCUMENT

 */
{

  //computes residues in arcsec
  residue  = (*rtc.slopes_res)(slrange(rtc.its),);
  residue -= residue(,avg);

  /////////////////////////////////
  // .... Computing covariance matrices
  /////////////////////////////////////////////////////
    
  //TS measurements covariance matrix
  Cee  = residue(,+) * residue(,+)/dimsof(residue)(0);

  // noise matrix
  Cnn  = *covMatrix.noiseCL;
  Cnn  = Cnn(slrange(rtc.its),slrange(rtc.its));

  // aliasing matrix
  Caa = *covMatrix.aliasing;
  
  if(rtc.obsMode == "MOAO"){
    R    = *rtc.R;
    Crr  = Caa(slrange(rtc.its),slrange(rtc.its));
    Cba  = Caa(slrange(rtc.its),norange(rtc.its));
    Cba  = R(,+)  * Cba(,+);
    Crr -= 2*Cba;
  }else{
    Crr  = -Caa(slrange(rtc.its),slrange(rtc.its));
  }
    
  //true residual covariance matrix estimation
  Cres = Cee - Cnn - Crr;
   
  /////////////////////////////////
  // .... Determining the residual OTF
  /////////////////////////////////////////////////////

  resOTF =  computeOTFtomographic(averageMode,Cee=Cres,verb=verb);

  //Final OTF
  resOTF *= roll(*otf.fit) * roll(*otf.static);

  return resOTF;
}



func displayPracResults(psf,otf,tel,atm,cam,budget)
/* DOCUMENT displayPracResults,psf,otf,tel,atm,cam,budget

   Displays results from PSF reconstruction performed by PRAC.
*/
{

  extern reDisp;
  
  if(reDisp){
    winclose;
    window,0;
    window,1;
    window,2;
    window,3;
    window,4,style="aanda.gs",dpi=90;
    window,5,style="aanda.gs",dpi=90;
    window,6,style="aanda.gs",dpi=90;
    window,7,style="aanda.gs",dpi=90;
  }


  meth = "Estimated";
  if(psfr == "ts")
    meth = "TS-based";
  else if(psfr == "rtc")
    meth = "RTC-based";

  
  nc = 64;
  nm = (cam.nPixelsCropped - nc)/2+1;
  np = (cam.nPixelsCropped + nc)/2;
  
  l = cam.nPixelsCropped * cam.pixSize;
  f = abs;


  // ....................... PSFS ...................................//
  
  window,0; clr;logxy,0,0; pli, f(*psf.ncpa)(nm:np,nm:np),-l/2,-l/2,l/2,l/2;
  pltitle,"Best bench PSF\n Strehl = " + var2str(arrondi(100*psf.SR_ncpa,1))+"%";
  xytitles,"Arcsec","Arcsec";
  pdf,"results/" + timedata + "_ncpaPSF.pdf";
    
  window,1; clr;pli, f(*psf.res)(nm:np,nm:np),-l/2,-l/2,l/2,l/2,cmax = f(psf.SR_sky);
  pltitle,meth +" PSF\n Strelh = " + var2str(arrondi(100*psf.SR_res,1))+"%";
  xytitles,"Arcsec","Arcsec";
  pdf,"results/" + timedata + "_" + psfr + "PSF" + ".pdf";
    
  window,2; clr;pli, f(*psf.sky)(nm:np,nm:np),-l/2,-l/2,l/2,l/2,cmax = f(psf.SR_sky);
  pltitle,"On-sky PSF\n Strehl = " + var2str(arrondi(100*psf.SR_sky,1)) +"%";
  xytitles,"Arcsec","Arcsec";
  //colorBar,min(log(abs(*psf.sky))),log(psf.SR_sky);
  colorBar,0,arrondi(100*psf.SR_sky,1);
  pdf,"results/" + timedata + "_skyPSF.pdf";
    
  window,3; clr; pli, f(*psf.diff)(nm:np,nm:np),-l/2,-l/2,l/2,l/2,cmax = f(.1*psf.SR_sky);
  pltitle,"Residue";
  xytitles,"Arcsec","Arcsec";
  colorBar,0,arrondi(10*psf.SR_sky,1);
  pdf,"results/" + timedata + "_" + psfr + "resPSF" + ".pdf";

  // .................... ERROR BREAKDOWN .........................//
    
  window,4; clr;
  y = [budget.res,budget.tomo,budget.alias,budget.noise,budget.bw,
       budget.fit,budget.ol,budget.static,budget.ncpa];
  labs = ["!s_!e","!s_Tomography","!s_Alias","!s_Noise",
          "!s_Servo", "!s_Fitting","!s_Go-to","!s_Static","!s_Ncpa"];

  plotsBarDiagram,y,labs,col1=[char(241)],title=1,step=2,thick=.1;
  pdf,"results/" + timedata + "_budget" + ".pdf";

  // .................... ENSQUARED ENERGY .........................//
    
  window,5; clr;gridxy,1,1;
  Edif = 100*sqrt(abs(*psf.EE_res^2 - *psf.EE_sky^2))/(*psf.EE_sky);
  m = max(Edif);
  plg, Edif, boxsize,marks=0;
  fcut = radian2arcsec*cam.lambda/tel.pitch;
  plg, [0,1.05*m],[fcut,fcut],type=2,marks=0;
  plg, [0,1.05*m],[1,1]*1.22*radian2arcsec*cam.lambda/tel.diam,type=2,marks=0;
  xytitles,"Integrated box width [arcsec]","Relative error on EE [%]";
  plt,"DM cut-off frequency",fcut*1.1,.8*m,tosys=1;
  plt,"1.22 !l/D",1.22*radian2arcsec*cam.lambda/tel.diam*1.05,10,tosys=1;
  range,0,1.05*m;
  limits,-0.1,max(boxsize)*1.05;
  pdf,"results/" + timedata + "_" + psfr + "EE" + ".pdf";

  // .................... RADIAL OTF AVERAGE .........................//
  /*
    window,6; clr;logxy,0,1;gridxy,1,1;
    os = circularAveragePsf(*otf.sky);os/=max(os);
    or = circularAveragePsf(*otf.res);or/=max(or);
    //defining the telescope OTF at the low resolution
    n    = cam.nPixelsCropped;
    otel = roll(pracTelescopeOTF(tel.obs,n));
    ot   = (otel/max(otel))(n/2+1,n/2+1:);
    dl   =  span(0,n/2,n/2) * tel.foV/tel.nPixels;

    plg,os,dl;
    plg,or,dl;
    plg,ot,dl,type=2,marks=0;
    
    xytitles,"Normalized frequency [D/!l units]","Circularly averaged OTF";
    plt,"Dashed line: Perfect telescope",.1,1e-2,tosys=1;
    plt,"A: Sky OTF",.1,.25e-2,tosys=1;
    plt,"B: "+meth+" OTF",.1,1e-3,tosys=1;
    range,1e-4,1.1;limits,-.1,1.1;
    pdf,"results/" + timedata + "_" + psfr + "OTF" + ".pdf";
  */
  // .................... TURBULENT PROFILES .........................//
    
  if(strpart(rtc.aoMode,1:4) == "MOAO"){
    window,7; clr;
    displayLayers,atm.cnh,atm.altitude,col=[char(242)],percent=1,thick=.4;
    displayLayers,-(*rtc.skyProfile)(,1),(*rtc.skyProfile)(,2),col=[char(241)],percent=1;
    limits,-100,100;
    range,-1,20;
    plt,"Calibration on-sky profile",-80,15,tosys=1;
    plt,"Post retrieved profile",30,15,tosys=1,color=[128,128,128];
    pdf,"results/" + timedata + "_calib" + ".pdf";
  }

}


func printPracResults(psf,budget,t_init,t_budg,t_rec,t_proc,tf)
/* DOCUMENT
   
 */
{
  asky = psf.moffat_sky;
  ares = psf.moffat_res;
  
  write,"-------------------------------------------------------";
  write,"Reconstruction performance";
  write,"-------------------------------------------------------";
    
  write,format="Sky  SR               = %.3g%s\n", 100*psf.SR_sky,"%";
  write,format="Rec. SR               = %.3g%s\n", 100*psf.SR_res,"%";
  write,format="Sky  FWHM             = %.4g%s\n", 1e3*psf.FWHM_sky," mas";
  write,format="Rec. FWHM             = %.4g%s\n", 1e3*psf.FWHM_res," mas";
  write,format="Chi^2                 = %.5g\n", psf.chi2;
  write,format="(diff/sky)(avg)       = %.3g\n", psf.diffAvg;
  write,format="(diff/sky)(rms)       = %.3g\n", psf.diffRms;
  write,format="Sky  Moffat profile   = %.3g,%.3g,%.3g,%.3g,%.3g\n",asky(1),asky(2),asky(3),asky(4),asky(5);
  write,format="Rec. Moffat profile   = %.3g,%.3g,%.3g,%.3g,%.3g\n",ares(1),ares(2),ares(3),ares(4),ares(5);
    
  write,"-------------------------------------------------------";
  write,"Strehl ratios breakdown (Available in estimation method)";
  write,"-------------------------------------------------------";
    
  write,format="SR Mar. from budget   = %.3g%s\n", 100*budget.SRmar,"%";
  write,format="SR Par. from budget   = %.3g%s\n", 100*budget.SRpar,"%";
  write,format="SR Bor. from budget   = %.3g%s\n", 100*budget.SRborn,"%";
  write,format="SR fit                = %.3g%s\n", 100*psf.SR_fit,"%";
  write,format="SR bw                 = %.3g%s\n", 100*psf.SR_bw,"%";
  write,format="SR tomo+alias+noise   = %.3g%s\n", 100*psf.SR_tomo,"%";
  write,format="SR static             = %.3g%s\n", 100*psf.SR_stats,"%";
  write,format="SR ncpa               = %.3g%s\n", 100*psf.SR_ncpa,"%";

  write,"-------------------------------------------------------";
  write,"WF error breakdown";
  write,"-------------------------------------------------------";
    
  write,format="Residual error        = %.4g nm rms\n", budget.res;
  write,format="Tomographic error     = %.4g nm rms\n", budget.tomo;
  write,format="Aliasing error        = %.4g nm rms\n", budget.alias;
  write,format="Noise error           = %.4g nm rms\n", budget.noise;
  write,format="Bandwidth error       = %.4g nm rms\n", budget.bw;
  write,format="Fitting error         = %.4g nm rms\n", budget.fit;
  write,format="Go-to error           = %.4g nm rms\n", budget.ol;
  write,format="Static error          = %.4g nm rms\n", budget.static;
  write,format="NCPA error            = %.4g nm rms\n", budget.ncpa;

  write,"-------------------------------------------------------";
  write,"Computation time";
  write,"-------------------------------------------------------";

  write,format="Instantiation      done on %.3g s\n",t_init;
  write,format="Error breakdown    done on %.3g s\n",t_budg;
  write,format="PSF reconstruction done on %.3g s\n",t_rec;
  write,format="PSF processing     done on %.3g s\n",t_proc;
  write,format="Total              done on %.3g s\n ",tf;
  
  reDisp = 0;
}



func concatenatePracResults(method,writeRes=)
/* DOCUMENT res = concatenatePracResults(method,writeRes=)

   Returns a pointer concatenating all results from PRAC using
   reconstruction method given by input "method".

 */
{
  pracResults = array(pointer,12);

  /////////////////////
  // ..... DATA IDENTITY + PARAMETERS IDENTIFICATION
  //////////////////////////////////////////////////////////
  
  //Data identity
  pracResults(1) = &strchar([timedata,rtc.aoMode,rtc.obsMode,rtc.recType]);
  //global parameters
  if(is_void(tf)) tf = 0;
  pracResults(2) = &[atm.r0,atm.L0,atm.v,tf];
  //turbulence profile
  pracResults(3) = &[atm.cnh,atm.altitude,atm.l0h,atm.vh];
  //tracking
  pracResults(4) = &[sys.tracking];
  //system parameters
  pracResults(5) = &[wfs.x,wfs.y,sys.xshift,sys.yshift,sys.magnification,sys.theta,sys.centroidGain];


  /////////////////////
  // ..... IMAGES
  //////////////////////////////////////////////////////////

  if(!is_void(*psf.sky)){
    //PSFs
    pracResults(6) = &[*psf.sky,*psf.res,*psf.diff];
    //Ensquared Energy
    pracResults(7) = &[*psf.EE_sky,*psf.EE_res];
    //OTFs
    pracResults(8) = &[*otf.sky,*otf.res];
  }
  
  /////////////////////
  // ..... PERFORMANCE
  //////////////////////////////////////////////////////////
  
  //error budget
  pracResults(9) = &[budget.res,budget.fit,budget.bw,budget.tomo,budget.noise,budget.alias,budget.static,budget.ncpa,budget.ol];
  //VED
  pracResults(10) = &budget.ved;
  //Strehl ratios
  pracResults(11) = &([psf.SR_sky,psf.SR_res,budget.SRmar,budget.SRpar,budget.SRborn,psf.SR_tomo,psf.SR_fit,psf.SR_bw,psf.SR_stats,psf.SR_ncpa,psf.FWHM_sky,psf.FWHM_res,psf.chi2]);
  //Moffat profiles
  pracResults(12) = &[psf.moffat_sky,psf.moffat_res];

  if(writeRes){
    savingDir = "results/pracResults/resultsWith" + method + "_Method/";
    if(!direxist(savingDir))
      system,"mkdir " + savingDir;
    writefits, savingDir + "pracResults_" + method + "_" + strpart(procDir,1:10) + "_" + timedata + ".fits",pracResults;

  }
  
  return pracResults;
}



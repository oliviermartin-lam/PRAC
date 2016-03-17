include, "pracDm.i"; 
include, "pracBwfit.i"; 
include, "pracTomography.i";
include, "pracLmfit.i";
include, "pracLearn.i";
include, "pracErrorbreakdown.i";
include, "pracAtmosphere.i";
include, "pracDirectPsfr.i";
include, "pracMakeOTF.i";

include, "canaryUtils.i";
include, "telescopeUtils.i"; 
include, "noiseUtils.i";
include, "imageUtils.i";
include, "contextUtils.i";
include, "fitsUtils.i";
include, "mathUtils.i";
include, "displayUtils.i";
include, "constantUtils.i";
include, "zernikeUtils.i";

if(!runDone){
  winclose;
  window,0;
  window,1;
  window,2;
  window,3;
  window,4,style="aanda.gs",dpi=90;
  window,5,style="aanda.gs",dpi=90;
  window,6,style="aanda.gs",dpi=90;
  window,7,style="aanda.gs",dpi=90;
  window,8,style="aanda.gs",dpi=90;
 }

func pracMain(timedata,Dir=,psfrMethod=,averageMode=,verb=,disp=,writeRes=)
/* DOCUMENT res = pracMain(timedata,Dir=,psfrMethod=,averageMode=,verb=,disp=,writeRes=)

  Launches the PSF reconstruction for processing the CANARY on-sky data slopestl acquired at "timedata"
  and storaged into the directory specified by the keyword "Dir"
  using the method given by the keyword "psfrMethod."

  INPUTS:

  - timedata (string)   : local hour of the slopestl file
  - Dir (string)        : Night directory in which the slopestl are storaged.
  It depends on how you've nammed your own data directories.
  - psfrMethod (string) : Reconstruction method: "estimation" for analytic reconstruction,
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
  
  res = pracMain("00h15m36s",Dir="2013_09_13_onsky/",verb=1,disp=1,psfrMethod = "estimation")

  performs the PSF-R using the estimation method on the slopestl file acquired at 00h15m36s 
  and storaged into the night directory "2013_09_13_onsky/".

  SEE ALSO: concatenatePracResults()
 */
{

  extern runDone;

  if(psfrMethod != "estimation" & psfrMethod != "ts" & psfrMethod != "rtc"){
    write,"Please make your choice between estimation, TS and RTC PSF-R methods";
    return [];
  }
  tic,10;
  include, "pracConfig.i",1;
  if(Dir){
    procDir = Dir;
    dataDirRoot = dataDir + procDir;
  }

  geometry = "square";
  if(is_void(psfrMethod))
    psfrMethod = "estimation";
  if(is_void(averageMode))
    averageMode = "Vii";

  
  /////////////////////
  //.... Initializing the data struct
  //////////////////////////////////////////////////////////
  
  include, "pracStructConfig.i",1;
  define_structs,timedata,verb=verb;

  /////////////////////
  // .... perfect telescope OTF
  //////////////////////////////////////////////////////////
  
  OTF_tel = OTF_telescope(tel.diam,tel.obs,tel.nPixels,tel.pixSize);
  otf.tel = &roll(OTF_tel);
  
  
  /////////////////////
  // .... OTF from fitting error
  //////////////////////////////////////////////////////////
   
  //Fitting OTF
  OTF_fit       = computeOTFfitting(geometry,verb=verb);
  otf.fit       = &roll(OTF_fit);
  psf.SR_fit    = sum((*otf.fit) * (*otf.tel));

  /////////////////////
  // .... OTF from static aberrations
  //////////////////////////////////////////////////////////

  //telescope, NCPA included
  OTF_telstats  = computeOTFstatic(PSF_stats,PSF_ncpa_fit,SR_stats,SR_ncpa);

  psf.SR_stats  = SR_stats;
  otf.static    = &roll(OTF_telstats);

  //Getting the best bench PSF
  OTF_ncpa      = getOTFncpa(cam.nPixelsCropped,procDir,SR_bench,PSF_ncpa,disp=disp);
  psf.ncpa      = &PSF_ncpa;
  psf.SR_ncpa   = budget.SRncpa = SR_ncpa;
  budget.ncpa   = 144;//sr2var(SR_ncpa,cam.lambda);

  /////////////////////
  // .... Tip-tilt OTF
  //////////////////////////////////////////////////////////

  OTF_tilt      = computeOTFtiptilt([],[0,0,0],Dtt);
  otf.tilt      = &OTF_tilt;
  psf.SR_tilt   = sum((*otf.tilt) * (*otf.tel));
  
 
  if(psfrMethod == "estimation"){
    
    /////////////////////
    // .... OTF from response delay time of the system
    //////////////////////////////////////////////////////////
  
    OTF_bw   = computeOTFbandwidth(geometry,rtc.obsMode,verb=verb);

    /////////////////////
    // .... OTF from tomographic residue, mode = "Uij", "Vii" or "intersample"
    //////////////////////////////////////////////////////////////////////////////
  
    OTF_tomo = computeOTFtomographic(averageMode,para=0,verb=verb);

       
    /////////////////////
    // .... OTF at the TS location
    //////////////////////////////////////////////////////////
  
    OTF_res =  OTF_telstats * OTF_fit * OTF_bw * OTF_tomo ;

       
    //telescope included into OTF_stats

    PSF_res = roll( fft(OTF_res).re );

    // .... Storage Strehl ratios
    psf.SR_tomo  = sum(OTF_tomo * OTF_tel);
    psf.SR_bw    = sum(OTF_bw * OTF_tel);
   
    // .... Storage OTFs
    otf.bw       = &roll(OTF_bw);
    otf.tomo     = &roll(OTF_tomo);
    otf.ncpa     = &roll(OTF_ncpa);
    
  
  }else if(psfrMethod == "ts"){

    //computes residues in arcsec
    residue  = (*rtc.slopes_res)(slrange(rtc.its),);
    //Split static and dynamic part
    stat     = residue(,avg);
    residue -= stat;

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
      Crr = computeCovSlopesError(Caa, *rtc.R);
      Caa(slrange(rtc.its),) = 0;
      Caa(,slrange(rtc.its)) = 0;
      Crr += computeCovSlopesError(Caa, *rtc.R);
    }else{
      Crr  = Caa(slrange(rtc.its),slrange(rtc.its));
    }
    //true residual covariance matrix estimation
    Cres = Cee - Cnn - Crr;
   
    /////////////////////////////////
    // .... Determining the residual OTF
    /////////////////////////////////////////////////////
    
    OTF_res =  computeOTFtomographic(averageMode,Dphi_res,Cee=Cres,verb=verb);

    //Final OTF
    OTF_res *= roll(*otf.fit) * roll(*otf.static) * roll(*otf.tilt);
    //grabbing the PSF
    PSF_res = roll(fft(OTF_res).re);

  }else if(psfrMethod == "rtc"){

    if(rtc.obsMode == "MOAO"){
      
      /////////////////////////////////
      // .... Computing the best MMSE tomographic reconstructor
      /////////////////////////////////////////////////////

      //grabbing matrices
      Coffoff  = (*covMatrix.parallel)(norange(rtc.its),norange(rtc.its));
      Conoff   = (*covMatrix.parallel)(slrange(rtc.its),norange(rtc.its));
      Conon    = (*covMatrix.parallel)(slrange(rtc.its),slrange(rtc.its));
      R        = *covMatrix.R;
            
      //Determination of the tomographic reconstruction residue covariance matrix
      tmp = R(,+) * Conoff(,+);
      Cdd = Conon - tmp - transpose(tmp) + R(,+)*(Coffoff(,+)*R(,+))(+,);
      
      /////////////////////////////////
      // .... Residual TS slopes from MMSE estimation
      /////////////////////////////////////////////////////

      // MMSE reconstruction
      soff = (*rtc.slopes_res)(norange(rtc.its),);
      son  = R(,+) * soff(+,); 

      // determining the delayed voltages
      V      = *rtc.volts;
      retard = 0.003*rtc.Fe + 1.05;
      fr     = int(retard);
      coef   = retard%1;
      Vconv  = coef*roll(V,[0,1+fr]) + (1.-coef)*roll(V,[0,fr]);
    
      sdm = (*rtc.mi)(,+) * Vconv(+,);
      sdm = mirror_SH7(sdm, wfs(rtc.its).sym);
      res = (son + sdm)(,3:-2);
      res-= res(,avg);

      //Residue covariance matrix
      Cres  = res(,+) * res(,+)/dimsof(res)(0) + Cdd;

      /////////////////////////////////
      // .... Determining the residual OTF
      /////////////////////////////////////////////////////
    
      OTF_res  = computeOTFtomographic(averageMode,Dphi_res,Cee=Cres,verb=verb);
      
    }else if(rtc.obsMode == "SCAO"){

           
      //Grabbing the voltages vectors in volts
      V    = (*rtc.volts)(,dif);
      //Command matrix in volts/arcsec
      MC   = *rtc.mc;      
      // noise matrix
      Cnn  = (*covMatrix.noiseCL)(slrange(rtc.its),slrange(rtc.its));
      Cnn  = MC(,+) * (Cnn(,+)*MC(,+))(+,);
      // aliasing matrix
      Caa  = (*covMatrix.aliasing)(slrange(rtc.its),slrange(rtc.its));
      Caa  = MC(,+) * (Caa(,+)*MC(,+))(+,);
      //Computing the volts covariance matrix
      Cvv  = V(,+) * V(,+)/dimsof(V)(0);
      Cvv  = Cvv - Cnn - Caa;
      
      //Computing the OTF
      OTF_res = makeOTF(Cvv,dphi,averageMode=averageMode,verb=verb);
     
    }
    
    //Final OTF
    OTF_res *= roll(*otf.fit) * roll(*otf.static) ;
    //grabbing the PSF
    PSF_res = roll(fft(OTF_res).re);   
  }
  
  /////////////////////
  // ..... Reconstructed PSF
  //////////////////////////////////////////////////////////

  nm = (tel.nPixels - cam.nPixelsCropped)/2+1;
  np = (tel.nPixels +  cam.nPixelsCropped)/2;

  psf.res    =  &PSF_res;//&PSF_res(nm:np,nm:np);
  *psf.res  /= sum(*psf.res);
  *psf.res  /= tel.airyPeak;
  psf.SR_res = max(*psf.res);
  psf.res    =  &(*psf.res)(nm:np,nm:np);
  otf.res    = &roll(fft(roll(*psf.res)).re);

  psf.ncpa_fit  = &(PSF_ncpa_fit(nm:np,nm:np));
   
  /////////////////////
  // .... Processing the sky PSF
  //////////////////////////////////////////////////////////

  //to set the maximum intensity at the middle pixel
  OTF_sky = roll(*otf.sky);
  
  PSF_sky  = roll(fft(OTF_sky).re);
  PSF_sky /= sum(PSF_sky);
  PSF_sky /= tel.airyPeak;
  
  
  psf.SR_sky   = max(PSF_sky);
  budget.SRsky = max(PSF_sky);

  psf.sky = &PSF_sky;
  
  
  // ... differential PSF
  
  psf.diff    = &(*psf.sky - *psf.res);
  w           = where(*psf.sky >0);
  psf.diffAvg = ((*psf.diff)(w)/ (*psf.sky)(w))(avg);
  psf.diffRms = ((*psf.diff)(w) /(*psf.sky)(w))(rms);
  psf.chi2    = sum(((*psf.diff)(w))^2/(*psf.sky)(w));
  
  /////////////////////
  // .... Error breakdown computation
  //////////////////////////////////////////////////////////
  
  PRAC_errorbreakdown,verb=verb;
  
  /////////////////////
  // ..... Getting Ensquared Energy and FWHM on both reconstructed/sky PSF
  //////////////////////////////////////////////////////////
  
  boxsize = EE = EE_sky = span(1,cam.nPixelsCropped-3.,cam.nPixelsCropped) * cam.pixSize;
  for(i=1; i<=cam.nPixelsCropped; i++){
    EE(i)     = getEE( 100*(*psf.res)/sum(*psf.res), cam.pixSize, boxsize(i));
    EE_sky(i) = getEE(100*(*psf.sky)/sum(*psf.sky), cam.pixSize, boxsize(i));
  }

  // Computes the EE at 10 lambda/D
  nee     = 10;
  dl      = indgen(cam.nPixelsCropped/2) * tel.foV/tel.nPixels;
  nn      = where(dl>= nee*radian2arcsec*cam.lambda/tel.diam)(1);
  EEsky_n = EE_sky(nn);
  EEres_n = EE(nn);

  // Computes the energy in the wings
  psf.skyWings = getWingsEnergy(*psf.sky,1.22,tel.pixSize,cam.lambda,tel.diam);
  psf.resWings = getWingsEnergy(*psf.res,1.22,tel.pixSize,cam.lambda,tel.diam);
    
  psf.EE_res   = &EE;
  psf.EE_sky   = &EE_sky;
  psf.FWHM_sky = getPsfFwhm(*psf.sky,tel.pixSize,asky,fit=2);
  psf.FWHM_res = getPsfFwhm(*psf.res,tel.pixSize,ares,fit=2);
  psf.skyMoffatProf = asky;
  psf.resMoffatProf = ares;
  
  /////////////////////
  // .... Display and verbose
  //////////////////////////////////////////////////////////
  
  if(disp){
    meth = "Estimated";
    if(psfrMethod == "ts")
      meth = "TS-based";
    else if(psfrMethod == "rtc")
      meth = "RTC-based";

    l = cam.nPixelsCropped * cam.pixSize;

    window,0; clr;logxy,0,0; pli, log(abs(*psf.ncpa)),-l/2,-l/2,l/2,l/2;
    pltitle,"Best bench PSF (log scale)\n Strehl = " + var2str(arrondi(100*psf.SR_ncpa,1))+"%";
    xytitles,"Arcsec","Arcsec";
    pdf,"results/" + timedata + "_ncpaPSF.pdf";
    
    window,1; clr;pli, log(abs(*psf.res)),-l/2,-l/2,l/2,l/2,cmax = log(psf.SR_sky);
    pltitle,meth +" PSF  (log scale)\n Strelh = " + var2str(arrondi(100*psf.SR_res,1))+"%";
    xytitles,"Arcsec","Arcsec";
    pdf,"results/" + timedata + "_" + psfrMethod + "PSF" + ".pdf";
    
    window,2; clr;pli, log(abs(*psf.sky)),-l/2,-l/2,l/2,l/2,cmax = log(psf.SR_sky);
    pltitle,"On-sky PSF  (log scale)\n Strehl = " + var2str(arrondi(100*psf.SR_sky,1)) +"%";
    xytitles,"Arcsec","Arcsec";
    colorBar,min(log(abs(*psf.sky))),log(psf.SR_sky);
    pdf,"results/" + timedata + "_skyPSF.pdf";
    
    window,3; clr; pli, log(abs(*psf.diff)),-l/2,-l/2,l/2,l/2,cmax = log(psf.SR_sky);
    pltitle,"Residue (log scale) ";
    xytitles,"Arcsec","Arcsec";
    pdf,"results/" + timedata + "_" + psfrMethod + "resPSF" + ".pdf";
    
    window,4; clr;
    y = [budget.res,
         budget.tomo,
         budget.alias,
         budget.noise,
         budget.bw,
         budget.fit,
         budget.ol,
         budget.static,
         budget.ncpa];
    labs = ["!s_!e",
            "!s_tomo",
            "!s_alias",
            "!s_noise",
            "!s_bw",
            "!s_fit",
            "!s_ol",
            "!s_static",
            "!s_ncpa"];

    plotsBarDiagram,y,labs,col1=[char(241)],title=1;
    pdf,"results/" + timedata + "_budget" + ".pdf";
    
    window,5; clr;
    plg, *psf.EE_res, boxsize,color=[128,128,128];
    plg, *psf.EE_sky, boxsize;
    plg, [100,100],[-0.1,max(boxsize)*1.05],type=2,marks=0;
    fcut = radian2arcsec*cam.lambda/tel.pitch;
    plg, [0,100],[fcut,fcut],type=2,marks=0;
    plg, [0,100],[1,1]*1.22*radian2arcsec*cam.lambda/tel.diam,type=2,marks=0;
    xytitles,"Angular separation from center (arcsec)","Ensquared Energy (%)";
    plt,"DM cut frequency",fcut*1.05,50,tosys=1;
    plt,"1.22 !l/D",1.22*radian2arcsec*cam.lambda/tel.diam*1.05,10,tosys=1;
    plt,"A: "+ meth + " PSF",1.05*fcut,30,tosys=1,color=[128,128,128];
    plt,"B: On-sky PSF",1.05*fcut,25,tosys=1;
    range,0,105;
    limits,-0.1,max(boxsize)*1.05;
    pdf,"results/" + timedata + "_" + psfrMethod + "EE" + ".pdf";
    
     window,6; clr;
    
    os = circularAveragePsf(*otf.sky);os/=max(os);
    or = circularAveragePsf(*otf.res);or/=max(or);
    //defining the telescope OTF at the low resolution
    otel = roll(OTF_telescope(tel.diam,tel.obs,cam.nPixelsCropped,\
                              tel.pixSize*tel.nPixels/cam.nPixelsCropped));
    n = dimsof(otel)(0);
    ot = (otel/max(otel))(n/2+1,n/2+1:);
    n = numberof(os);
    dl =  span(0,cam.nPixelsCropped/2,n) * tel.foV/tel.nPixels;
    plg,os,dl;
    plg,or,dl;
    plg,ot,dl,type=2,marks=0;
    xytitles,"!r/D","Normalized radially averaged OTF";
    plt,"Dashed line: Perfect telescope",1,1,tosys=1;
    plt,"A: Sky OTF",1,0.8,tosys=1;
    plt,"B: "+meth+" OTF",1,0.6,tosys=1;
    range,-.1,1.1;
    pdf,"results/" + timedata + "_" + psfrMethod + "OTF" + ".pdf";

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

    window,8; clr;
    clr;
    y  = atm.altitude/1000.;
    x1 = atm.cnh;
    x2 = atm.l0h;
    x3 = atm.vh;

    //managing outer scale
    w = where(x2>100);
    if(is_array(w))
      x2(w) = 0.
        
    nn = sort(y);
    y  = y(nn);
    x1 = x1(nn);
    x2 = x2(nn);
    x3 = x3(nn);
    thick = .15;
    dy = 2*thick;

    for(i=1; i<=atm.nLayers; i++) {
      plfp,[char(242)], [1,1,1,1]*y(i)+dy+[-1,-1,1,1]*thick,[0,1,1,0]*x1(i), [4],edges=0;
      plfp,[char(243)], [1,1,1,1]*y(i)+[-1,-1,1,1]*thick,[0,1,1,0]*x2(i), [4],edges=0;
      plfp,[char(241)], [1,1,1,1]*y(i)-dy+[-1,-1,1,1]*thick,[0,1,1,0]*x3(i), [4],edges=0;
    }
    range,-1,20;
    m = max(grow(x1,x2,x3));
    limits,-1,m*1.1;
    plt,"r_0^-5/3^(h) (m^-5/3^)",m/2,19,color = [128,128,128],tosys=1;
    plt,"L_0_(h) (m)",m/2,18,color = [100,100,100],tosys=1;
    plt,"v(h) (m/s)",m/2,17,color = "black",tosys=1;
    xytitles,"Profiles","Altitude (km)";
    pdf,"results/" + timedata + "_profiles.pdf";
  }

  tf = tac(10);
  
  if(verb){

    write,"-------------------------------------------------------";
    write,"Reconstruction performance";
    write,"-------------------------------------------------------";
    
    write,format="Sky  SR               = %.3g%s\n", 100*psf.SR_sky,"%";
    write,format="Rec. SR               = %.3g%s\n", 100*psf.SR_res,"%";
    write,format="Sky  FWHM             = %.4g%s\n", 1e3*psf.FWHM_sky," mas";
    write,format="Rec. FWHM             = %.4g%s\n", 1e3*psf.FWHM_res," mas";
    write,format="Sky  EE at %dlambda/D = %.4g%s\n", nee,EEsky_n," %";
    write,format="Rec. EE at %dlambda/D = %.4g%s\n", nee,EEres_n," %";
    write,format="Sky  wings energy     = %.4g%s\n", psf.skyWings," %";
    write,format="Rec. wings energy     = %.4g%s\n", psf.resWings," %";
    write,format="Chi^2                 = %.3g\n", psf.chi2;
    write,format="(diff/sky)(avg)       = %.3g\n", psf.diffAvg;
    write,format="(diff/sky)(rms)       = %.3g\n", psf.diffRms;
    write,format="Sky  Moffat profile   = %.3g,%.3g,%.3g\n",asky(1),asky(2),asky(3);
    write,format="Rec. Moffat profile   = %.3g,%.3g,%.3g\n",ares(1),ares(2),ares(3);
    
    write,"-------------------------------------------------------";
    write,"Strehl ratios breakdown (Avaliable in estimation method)";
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
    write,format="PSF reconstruction done on %.3g s\n ",tf;
  }

  pracResults = concatenatePracResults(psfrMethod + "_" + averageMode,writeRes=writeRes);
  runDone = 1;
  return pracResults;
}


func concatenatePracResults(method,writeRes=)
/* DOCUMENT

 */
{
  pracResults = array(pointer,12);

  /////////////////////
  // ..... DATA IDENTITY + PARAMETERS IDENTIFICATION
  //////////////////////////////////////////////////////////
  
  //Data identity
  pracResults(1) = &strchar([timedata,rtc.aoMode,rtc.obsMode,rtc.recType]);
  //global parameters
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

  //PSFs
  pracResults(6) = &[*psf.sky,*psf.res,*psf.diff,*psf.ncpa];
  //Ensquared Energy
  pracResults(7) = &[*psf.EE_sky,*psf.EE_res];

  pracResults(8) = &[*otf.sky,*otf.res];
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
  pracResults(12) = &[psf.skyMoffatProf,psf.resMoffatProf];

  if(writeRes){
    savingDir = "results/pracResults/resultsWith" + method + "_Method/";
    if(!direxist(savingDir))
      system,"mkdir " + savingDir;
    writefits, savingDir + "pracResults_" + method + "_" + strpart(procDir,1:10) + "_" + timedata + ".fits",pracResults;

  }
  
  return pracResults;
}



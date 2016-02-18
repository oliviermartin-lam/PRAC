include, "pracDm.i"; 
include, "pracBwfit.i"; 
include, "pracTomography.i";
include, "pracLmfit.i";
include, "pracLearn.i";
include, "pracErrorbreakdown.i";
include, "pracAtmosphere.i";
include, "pracDirectPsfr.i";

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



func pracMain(timedata,psfrMethod=,averageMode=,Dir=,verb=,disp=,writeRes=)
/*DOCUMENT p = pracMain("00h15m36s",Dir="2013_09_13_onsky/",verb=1,disp=1,psfrMethod = "analytic");
  SCAO: "00h14m37s" GLAO:"00h16m30s"
  
 
 */
{

  tic,10;
  include, "pracConfig.i",1;
  if(Dir){
    procDir = Dir;
    dataDirRoot = dataDir + procDir;
  }

  geometry = "square";
  if(is_void(psfrMethod))
    psfrMethod = "analytic";
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
  //defining the telescope OTF at the low resolution
  otel    = roll(OTF_telescope(tel.diam,tel.obs,cam.nPixelsCropped,tel.pixSize*tel.nPixels/cam.nPixelsCropped));
  
  /////////////////////
  // .... OTF from fitting error + ncpa
  //////////////////////////////////////////////////////////

  //Getting the best bench PSF
  OTF_ncpa      = getOTFncpa(cam.nPixelsCropped,procDir,SR_ncpa,PSF_ncpa,disp=disp);
  budget.SRncpa = psf.SR_ncpa = SR_ncpa;
  psf.ncpa      = &PSF_ncpa;
  budget.ncpa   = sr2var(SR_ncpa,cam.lambda);
  
  //Fitting OTF
  power = 3;
  OTF_fit       = computeOTFfitting(geometry,SR_ncpa,power,verb=verb);//ncpa included
  otf.fit       = &roll(OTF_fit);
  psf.SR_fit    = sum((*otf.fit) * (*otf.tel))/sum(*otf.tel);

  /////////////////////
  // .... OTF from static aberration except NCPA
  //////////////////////////////////////////////////////////
  
  OTF_telstats  = computeOTFstatic(PSF_stats);//telescope included
  otf.static    = &roll(OTF_telstats);

  if(psfrMethod == "rtc" && rtc.obsMode == "SCAO"){
    psfrMethod = "ls";
    write,"Switch to ls method in SCAO case";
  }

  if(psfrMethod == "analytic"){
    
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
  
    OTF_res = OTF_fit * OTF_bw * OTF_tomo ;

    if(rtc.obsMode == "MOAO"){
      OTF_res *=  OTF_telstats;
    }else{
      OTF_res *= OTF_tel;
    }
    //telescope included into OTF_stats

    PSF_res = roll( fft(OTF_res).re );

    // .... Storage Strehl ratios
    psf.SR_tomo  = sum(OTF_tomo * OTF_tel)/sum(OTF_tel);
    psf.SR_bw    = sum(OTF_bw * OTF_tel)/sum(OTF_tel);
    psf.SR_stats = sum(OTF_telstats)/sum(OTF_tel);

    // .... Storage OTFs
    otf.bw       = &roll(OTF_bw);
    otf.tomo     = &roll(OTF_tomo);
    otf.ncpa     = &roll(OTF_ncpa);
    
  
  }else if(psfrMethod == "ls"){

    //computes residues in arcsec
    residue  = (*rtc.slopes_res)(slrange(rtc.its),);
    stat     = residue(,avg);
    residue -= stat;

    /////////////////////////////////
    // .... Computing covariance matrices
    /////////////////////////////////////////////////////
    
    //TS measurements covariance matrix
    Cee  = residue(,+) * residue(,+)/dimsof(residue)(0);

    // noise matrix
    Cnn  = *covMatrix.noise;
    Cnn  = Cnn(slrange(rtc.its),slrange(rtc.its));

    // aliasing matrix
    Caa = *covMatrix.aliasing;
    if(rtc.obsMode == "MOAO"){
      Crr = computeCovSlopesError(Caa,*rtc.R);
      Caa(slrange(rtc.its),) = 0;
      Caa(,slrange(rtc.its)) = 0;
      Crr += computeCovSlopesError(Caa,*rtc.R);
    }else{
      Crr = Caa(slrange(rtc.its),slrange(rtc.its));
      fact = 1.+ filteringNoiseFactor(rtc.loopGain,rtc.frameDelay,rtc.obsMode);
      Cnn *= fact;
    }
    //true residual covariance matrix estimation
    Cres = Cee - Cnn - Crr;
   
    /////////////////////////////////
    // .... Determining the residual OTF
    /////////////////////////////////////////////////////
    
    OTF_res =  computeOTFtomographic(averageMode,Dphi_res,Cee=Cres,verb=verb);

    //Final OTF
    OTF_res *= roll(*otf.fit) * roll(*otf.static);
    //grabbing the PSF
    PSF_res = roll(fft(OTF_res).re);

  }else if(psfrMethod == "rtc"){

    /////////////////////////////////
    // .... Computing the best MMSE tomographic reconstructor
    /////////////////////////////////////////////////////

    //grabbing matrices
    covLearn = *covMatrix.learn;
    covNoise = *covMatrix.noise;
    covPara = *covMatrix.parallel;
    covMeas = *covMatrix.slopes;
    covAlias = *covMatrix.aliasing;
    
    Coffoff = (covLearn +  covNoise)(norange(rtc.its),norange(rtc.its));
    Conoff  = covPara(slrange(rtc.its),norange(rtc.its));
    Conon   = covPara(slrange(rtc.its),slrange(rtc.its));

    //inversion
    iCoffoff = invgen(Coffoff,10);

    //computation
    R = Conoff(,+) * iCoffoff(+,);
        
    //Determination of the tomographic reconstruction residue covariance matrix
    tmp = R(,+) * Conoff(,+);
    Cdd = Conon - tmp - transpose(tmp) + R(,+)*(Coffoff(,+)*R(,+))(+,);
    
    //Unbiases from noise/aliasing propagated through reconstructor
    C   = (covNoise + covAlias)(norange(rtc.its),norange(rtc.its));
    Cnn =  R(,+)*(C(,+)*R(,+))(+,);
    
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
    Vconv  = coef*roll(V,[0,fr]) + (1.-coef)*roll(V,[0,fr+1]);
    
    sdm = (*rtc.mi)(,+) * Vconv(+,);
    sdm = mirror_SH7(sdm, wfs(rtc.its).sym);
    res = (son + sdm)(,3:-2);
    res-= res(,avg);

    
    /////////////////////////////////
    // .... Determining the residual OTF
    /////////////////////////////////////////////////////

    //Residue covariance matrix
    Cres  = res(,+) * res(,+)/dimsof(res)(0) + Cdd - Cnn;
    
    OTF_res =  computeOTFtomographic(averageMode,Dphi_res,Cee=Cres,verb=verb);

    //Final OTF
    OTF_res *= roll(*otf.fit) * roll(*otf.static);
    //grabbing the PSF
    PSF_res = roll(fft(OTF_res).re);

    
  }
  
  /////////////////////
  // ..... Reconstructed PSF
  //////////////////////////////////////////////////////////

  
  nm = (tel.nPixels - cam.nPixelsCropped)/2+1;
  np = (tel.nPixels +  cam.nPixelsCropped)/2;
    
  psf.res   =  &PSF_res(nm:np,nm:np);
  *psf.res  /= sum(*psf.res);
  *psf.res  /= tel.airyPeak;
  psf.SR_res = max(*psf.res);
  otf.res    = &roll(fft(roll(*psf.res)).re);

   
  /////////////////////
  // .... Processing the sky PSF
  //////////////////////////////////////////////////////////

  //to set the maximum intensity at the middle pixel
  OTF_sky = roll(*otf.sky);
  
  PSF_sky = roll(fft(OTF_sky).re);
  PSF_sky/= sum(PSF_sky);
  PSF_sky/= tel.airyPeak;
  
  psf.SR_sky   = max(PSF_sky);
  budget.SRsky = max(PSF_sky);

  psf.sky = &PSF_sky;
  
  
  // ... differential PSF
  
  psf.diff    = &(*psf.sky - *psf.res);
  psf.diffSum = sum(abs(*psf.diff));
  psf.diffRms = (*psf.diff)(rms);

  
  /////////////////////
  // .... Error breakdown computation
  //////////////////////////////////////////////////////////
  
  PRAC_errorbreakdown,verb=verb;
  
  /////////////////////
  // ..... Getting Ensquared Energy and FWHM on both reconstructed/sky PSF
  //////////////////////////////////////////////////////////
  
  boxsize = EE = EE_sky = span(1,cam.nPixelsCropped-3.,cam.nPixelsCropped) * cam.pixSize;
  for(i=1; i<=cam.nPixelsCropped; i++){
    EE(i) = getEE( 100*(*psf.res)/sum(*psf.res), cam.pixSize, boxsize(i));
    EE_sky(i) = getEE(100*(*psf.sky)/sum(*psf.sky), cam.pixSize, boxsize(i));
  }

  psfsky  = (*psf.sky)(cam.nPixelsCropped/2+1:,cam.nPixelsCropped/2+1);
  fwhmSky = getFWHM(psfsky/max(psfsky),1./tel.pixSize);
  psfres  = (*psf.res)(cam.nPixelsCropped/2+1:,cam.nPixelsCropped/2+1);
  fwhmRes = getFWHM(psfres/max(psfres),1./tel.pixSize);

  
  psf.EE_res   = &EE;
  psf.EE_sky   = &EE_sky;
  psf.FWHM_sky = fwhmSky;
  psf.FWHM_res = fwhmRes;

  
  /////////////////////
  // .... Display and verbose
  //////////////////////////////////////////////////////////
  
  if(disp){
    meth = "Analytic";
    if(psfrMethod == "ls")
      meth = "TS-based LS";
    else if(psfrMethod == "rtc")
      meth = "RTC-based reconstruction";

    l = cam.nPixelsCropped * cam.pixSize;

    
    window,0; clr;logxy,0,0; pli, *psf.ncpa,-l/2,-l/2,l/2,l/2,cmin=0;
    pltitle,"Best bench PSF with SR = " + var2str(arrondi(100*psf.SR_ncpa,1))+"%";
    xytitles,"Arcsec","Arcsec";
    pdf,"results/" + timedata + "_ncpaPSF.pdf";
    
    window,1; clr;pli, *psf.res,-l/2,-l/2,l/2,l/2,cmin=0,cmax = psf.SR_sky;
    pltitle,meth +" PSF with SR = " + var2str(arrondi(100*psf.SR_res,1))+"%";
    xytitles,"Arcsec","Arcsec";
    pdf,"results/" + timedata + "_" + psfrMethod + "PSF" + ".pdf";
    
    window,2; clr;pli, *psf.sky,-l/2,-l/2,l/2,l/2,cmin=0,cmax = psf.SR_sky;
    pltitle,"On-sky PSF with SR = " + var2str(arrondi(100*psf.SR_sky,1)) +"%";
    xytitles,"Arcsec","Arcsec";
    pdf,"results/" + timedata + "_skyPSF.pdf";
    
    window,3; clr; pli, *psf.diff,-l/2,-l/2,l/2,l/2,cmin=0;
    pltitle,"Residual of the reconstruction ";
    xytitles,"Arcsec","Arcsec";
    pdf,"results/" + timedata + "_" + psfrMethod + "resPSF" + ".pdf";
    
    winkill,4;window,4,style="aanda.gs",dpi=90;clr;
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
    
    winkill,5;window,5,style="aanda.gs",dpi=90;clr;
    plg, *psf.EE_res, boxsize,color=[128,128,128];
    plg, *psf.EE_sky, boxsize;
    plg, [100,100],[-0.1,max(boxsize)*1.05],type=2,marks=0;
    fcut = radian2arcsec*cam.lambda/tel.pitch;
    plg, [0,100],[fcut,fcut],type=2,marks=0;
    plg, [0,100],[1,1]*1.22*radian2arcsec*cam.lambda/tel.diam,type=2,marks=0;
    xytitles,"Angular separation from center (arcsec)","Ensquared Energy (%)";
    plt,"DM cut frequency",fcut*1.05,50,tosys=1;
    plt,"1.22 !l/D",1.22*radian2arcsec*cam.lambda/tel.diam*1.05,10,tosys=1;
    plt,"A: "+ meth + " PSF",0.5,30,tosys=1,color=[128,128,128];
    plt,"B: On-sky PSF",0.5,25,tosys=1;
    range,0,105;
    limits,-0.1,max(boxsize)*1.05;
    pdf,"results/" + timedata + "_" + psfrMethod + "EE" + ".pdf";
    
    winkill,6;window,6,style="aanda.gs",dpi=90;clr;
    dl =  indgen(cam.nPixelsCropped/2) * tel.foV/tel.nPixels;
    plg,(*otf.sky)(cam.nPixelsCropped/2+1,cam.nPixelsCropped/2+1:)/max(*otf.sky),dl;
    plg,(*otf.res)(cam.nPixelsCropped/2+1,cam.nPixelsCropped/2+1:)/max(*otf.res),dl;
    plg,otel(cam.nPixelsCropped/2+1,cam.nPixelsCropped/2+1:)/max(otel),dl,type=2,marks=0;
    xytitles,"D/!l","Normalized OTF";
    plt,"Dashed line: Perfect telescope",1,1,tosys=1;
    plt,"A: Sky OTF",1,0.8,tosys=1;
    plt,"B: "+meth+" OTF",1,0.6,tosys=1;
    range,-.1,1.1;
    pdf,"results/" + timedata + "_" + psfrMethod + "OTF" + ".pdf";

    if(strpart(rtc.aoMode,1:4) == "MOAO"){
      winkill,7;window,7,style="aanda.gs",dpi=90;clr;
      displayLayers,atm.cnh,atm.altitude,col=[char(242)],percent=1,thick=.4;
      displayLayers,-(*rtc.skyProfile)(,1),(*rtc.skyProfile)(,2),col=[char(241)],percent=1;
      limits,-100,100;
      range,-1,20;
      plt,"Calibration on-sky profile",-80,15,tosys=1;
      plt,"Post retrieved profile",30,15,tosys=1,color=[128,128,128];
      pdf,"results/" + timedata + "_calib" + ".pdf";
    }

    winkill,8;window,8,style="aanda.gs",dpi=90;
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
    write,format="SR on-sky           = %.3g%s\n", 100*psf.SR_sky,"%";
    write,format="SR reconstructed    = %.3g%s\n", 100*psf.SR_res,"%";
    write,format="SR Mar. from budget = %.3g%s\n", 100*budget.SRmar,"%";
    write,format="SR Par. from budget = %.3g%s\n", 100*budget.SRpar,"%";
    write,format="SR Bor. from budget = %.3g%s\n", 100*budget.SRborn,"%";
    write,format="SR fit              = %.3g%s\n", 100*psf.SR_fit,"%";
    write,format="SR bw               = %.3g%s\n", 100*psf.SR_bw,"%";
    write,format="SR tomo+alias+noise = %.3g%s\n", 100*psf.SR_tomo,"%";
    write,format="SR static           = %.3g%s\n", 100*psf.SR_stats,"%";
    write,format="SR ncpa             = %.3g%s\n", 100*psf.SR_ncpa,"%";
    write,format="Sum on recons.      = %.3g\n", psf.diffSum;
    write,format="Rms on recons.      = %.3g\n", psf.diffRms;
    write,format="Residual error      = %.4g nm rms\n", budget.res;
    write,format="Tomographic error   = %.4g nm rms\n", budget.tomo;
    write,format="Aliasing error      = %.4g nm rms\n", budget.alias;
    write,format="Noise error         = %.4g nm rms\n", budget.noise;
    write,format="Bandwidth error     = %.4g nm rms\n", budget.bw;
    write,format="Fitting error       = %.4g nm rms\n", budget.fit;
    write,format="Go-to error         = %.4g nm rms\n", budget.ol;
    write,format="Static error        = %.4g nm rms\n", budget.static;
    write,format="NCPA error          = %.4g nm rms\n", budget.ncpa;
    write,format="PSF reconstruction done on %.3g s\n",tf;
  }

  pracResults = concatenatePracResults(psfrMethod + "_" + averageMode,writeRes=writeRes);
  
  return pracResults;
}


func concatenatePracResults(method,writeRes=)
/* DOCUMENT

 */
{
  pracResults = array(pointer,11);

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
  pracResults(11) = &(100*[psf.SR_sky,psf.SR_res,budget.SRmar,budget.SRpar,budget.SRborn,psf.SR_tomo,psf.SR_fit,psf.SR_bw,psf.SR_stats,psf.SR_ncpa,psf.diffSum,psf.diffRms]);

  if(writeRes){
    savingDir = "results/pracResults/resultsWith" + method + "_Method/";
    if(!direxist(savingDir))
      system,"mkdir " + savingDir;
    writefits, savingDir + "pracResults_" + method + "_" + strpart(procDir,1:10) + "_" + timedata + ".fits",pracResults;

  }
  
  return pracResults;
}


/*

  else if(psfrMethod == "instantaneous"){

    // ................. TO BE DEEPLY REVIEWED ...............................//

    
    //computes residues in arcsec
    residue  = (*rtc.slopes_res)(slrange(rtc.its),);
    stats    = residue(,avg);
    dyn      = residue-stats;
    
    /////////////////////
    // .... Determining the denoised and dealiased residual slopes
    //////////////////////////////////////////////////////////

    gpara = mmseParallelModes(residue);

    /////////////////////
    // .... Average the Zernike modes every exposure time
    //////////////////////////////////////////////////////////

    zer2rad = pi*tel.diam/(radian2arcsec*cam.lambda); 
    //computing the Zernike modes in radians
    a_g  = (*sys.slopesToZernikeMatrix)(,+) * gpara(+,) * zer2rad;
    Zi   = readfits("fitsFiles/zernikeModes.fits",err=1);
    nzer = dimsof(Zi)(0);
    
    //constants
    nframes = dimsof(gpara)(0);
    t       = int(cam.exposureTime*rtc.Fe);
    nstep   = int(nframes/t);
    a_res_avg = a_stat_avg = array(0.,nzer,nstep+1);

    // Loop on exposure time step
    for(i = 1;i<=nstep+1;i++){
      //kth exposition
      if(i<=nstep)
        tk = 1 + (i-1)*t:i*t;
      else
        tk = nstep*t:0;
      //average on dynamic modes: the phase is averaged
      a_res_avg(,i)  = (a_g(,tk))(,avg);
      //static modes: the slopes are average to get the static phase
      stat_k = (residue(,tk))(,avg);
      a_stat_k = (*sys.slopesToZernikeMatrix)(,+) * stat_k(+);
      a_stat_avg(,i) = a_stat_k*zer2rad;
    }
    
    //Adding static

    a_res_avg  += a_stat_avg;

    /////////////////////
    // .... Determining the PSF
    ////////////////////////////////////////////

    P  = circularPupFunction(tel.diam,tel.obs,tel.nPixels,tel.pixSize);

    phi_res_k = PSF_res = 0*P;

    if(nstep != 1){
      for(k=1;k<=nstep+1;k++){
        phi_res_k *=0;
        for(i=1;i<=nzer;i++){
          phi_res_k += a_res_avg(i,k) * Zi(,,i);
        }
        //Electrical Field
        E = roll(P*exp(1i*phi_res_k));
        //PSF
        PSF_res += abs(fft(E))^2;
        write,format="Job done:%.3g%s\r",100.*k/(nstep+1.),"%";
      }
    }else{
      for(i=1;i<=nzer;i++){
        phi_res_k += a_res_avg(i,1) * Zi(,,i);
      }
      E = roll(P*exp(1i*phi_res_k));
      PSF_res = abs(fft(E))^2;
    }

    /////////////////////////////////
    // .... Adding high order modes and normalization
    /////////////////////////////////////////////////////
    
    // .... Adding fitting and computing Strehl ratio
    OTF_res = fft(PSF_res).re;
    OTF_res *= roll(*otf.fit);
    PSF_res = roll(fft(OTF_res).re);
  }
*/

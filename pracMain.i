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


pldefault,font="times";
pldefault,palette="earth.gp";
pltitle_font = "times";

func pracMain(timedata,psfrMethod=,averageMode=,Dir=,verb=,disp=,writeRes=)
/*DOCUMENT p = pracMain("02h31m34s",verb=1,disp=0); p = pracMain("02h30m48s",verb=1,disp=1);p = pracMain("02h31m11s",verb=1,disp=1,psfrMethod = "analytic",averageMode = "Vii",writeRes=1);
 
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
  otf.tel    = &roll(OTF_tel);
  otel = roll(OTF_telescope(tel.diam,tel.obs,cam.nPixelsCropped,tel.pixSize*tel.nPixels/cam.nPixelsCropped));
  
  /////////////////////
  // .... OTF from fitting error + ncpa
  //////////////////////////////////////////////////////////
  OTF_ncpa = getOTFncpa(cam.nPixelsCropped,procDir,SR_ncpa,PSF_ncpa,disp=disp);
  budget.SRncpa = psf.SR_ncpa = SR_ncpa;
  psf.ncpa     = &PSF_ncpa;
  
  OTF_fit = computeOTFfitting(geometry,SR_ncpa,verb=verb);//ncpa included
  otf.fit = &roll(OTF_fit);
  psf.SR_fit   = sum((*otf.fit) * (*otf.tel));

  /////////////////////
  // .... OTF from static aberration except NCPA
  //////////////////////////////////////////////////////////
  
  OTF_telstats  =  computeOTFstatic(PSF_stats);//telescope included
  otf.static    = &roll(OTF_telstats);

  if(psfrMethod == "analytic"){
    
    /////////////////////
    // .... OTF from response delay time of the system
    //////////////////////////////////////////////////////////
  
    OTF_bw = computeOTFbandwidth(geometry,rtc.obsMode,verb=verb);

    /////////////////////
    // .... OTF from tomographic residue, mode = "Uij", "Vii" or "intersample"
    //////////////////////////////////////////////////////////////////////////////
  
    OTF_tomo = computeOTFtomographic(averageMode,verb=verb);

    /////////////////////
    // .... OTF at the TS location
    //////////////////////////////////////////////////////////
  
    OTF_res = OTF_fit * OTF_bw * OTF_tomo * OTF_telstats ; //telescope included into OTF_stats

    PSF_res  = roll( fft(OTF_res).re );

    // .... Storage Strehl ratios
    psf.SR_tomo  = sum(OTF_tomo * OTF_tel);
    psf.SR_bw    = sum(OTF_bw * OTF_tel);
    psf.SR_stats = sum(OTF_telstats);

    // .... Storage OTFs
    otf.bw       = &roll(OTF_bw);
    otf.tomo     = &roll(OTF_tomo);
    otf.ncpa     = &roll(OTF_ncpa);
    
  
  }else if(psfrMethod == "mmse"){
   
    //computes residues in arcsec
    residue  = (*rtc.slopes_res)(slrange(rtc.its),);
    stats    = residue(,avg);

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
    Zi = readfits("fitsFiles/zernikeModes.fits",err=1);
    nzer = dimsof(Zi)(0);
    
    //constants
    t       = int(cam.exposureTime*rtc.Fe);
    nframes = dimsof(gpara)(0);
    nstep   = int(nframes/t);
    a_res_avg = a_stat_avg = array(0.,nzer,nstep+1);

    // Loop on exposure time step
    for(i = 1;i<=nstep;i++){
      tk = 1 + (i-1)*t:i*t;//kth exposition
      //average on dynamic modes: the phase is averaged
      a_res_avg(,i)  = (a_g(,tk))(,avg);
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
    a_res_avg(,0) =(a_g(,t*nstep:))(,avg);

    //Adding static 
    a_res_avg   += a_stat_avg;

    /////////////////////
    // .... Determining the PSF
    ////////////////////////////////////////////

    P  = circularPupFunction(tel.diam,tel.obs,tel.nPixels,tel.pixSize);

    phi_res_k = PSF_res = 0*P;
    
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

    /////////////////////////////////
    // .... Adding high order modes and normalization
    /////////////////////////////////////////////////////
    
    // .... Adding fitting + ncpa
    OTF_res = fft(PSF_res).re;
    OTF_res *= roll(*otf.fit);

    //grabbing the PSF
    PSF_res = roll(fft(OTF_res).re);

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
    Cnn  = 0*Cee;
    takesDiag(Cnn) = getNoiseVar(residue);
    // aliasing matrix
    Crr = *covMatrix.aliasing;
    if(rtc.obsMode == "MOAO")
      Crr = computeCovSlopesError(Crr,*rtc.R);
    else
      Crr = Crr(slrange(rtc.its),slrange(rtc.its));
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
  }
  
  /////////////////////
  // ..... Reconstructed PSF
  //////////////////////////////////////////////////////////

  //cropping
  nm = (tel.nPixels - cam.nPixelsCropped)/2+1;
  np = (tel.nPixels +  cam.nPixelsCropped)/2;
  psf.res =  &PSF_res(nm:np,nm:np);
  *psf.res /= sum(*psf.res);
  *psf.res /= tel.airyPeak;
  psf.SR_res = max(*psf.res);
  
  otf.res  = &roll(fft(roll(*psf.res)).re);
  
  /////////////////////
  // .... Processing the sky PSF
  //////////////////////////////////////////////////////////

  //to fixe the maximum Intensity at the middle
  OTF_sky = roll(*otf.sky);
  PSF_sky = roll(fft(OTF_sky).re);
  //normalization
  PSF_sky /= sum(PSF_sky);
  PSF_sky /= tel.airyPeak;
  psf.sky  = &PSF_sky;
  psf.SR_sky   = max(*psf.sky);
  
  // ... differential PSF
  
  psf.diff = &(abs(*psf.sky - *psf.res));
  psf.diffPsf = sum((*psf.diff)(*))/sum((*psf.sky)(*));

  
  /////////////////////
  // .... Error breakdown computation
  //////////////////////////////////////////////////////////
  
  PRAC_errorbreakdown,verb=verb;

  budget.ncpa  = sr2var(SR_ncpa,cam.lambda);
  budget.SRsky =  max(PSF_sky);
  
  
  /////////////////////
  // ..... Getting Ensquared Energy and FWHM on both reconstructed/sky PSF
  //////////////////////////////////////////////////////////
  
  boxsize = EE = EE_sky = span(1,cam.nPixelsCropped-3.,cam.nPixelsCropped) * cam.pixSize;
  for(i=1; i<=cam.nPixelsCropped; i++){
    EE(i) = getEE( 100*(*psf.res)/sum(*psf.res), cam.pixSize, boxsize(i));
    EE_sky(i) = getEE(100*(*psf.sky)/sum(*psf.sky), cam.pixSize, boxsize(i));
  }

  psfsky = (*psf.sky)(cam.nPixelsCropped/2+1:,cam.nPixelsCropped/2+1);
  fwhmSky = getFWHM(psfsky/max(psfsky),1./tel.pixSize);
  psfres = (*psf.res)(cam.nPixelsCropped/2+1:,cam.nPixelsCropped/2+1);
  fwhmRes = getFWHM(psfres/max(psfres),1./tel.pixSize);

  
  psf.EE_res   = &EE;
  psf.EE_sky   = &EE_sky;
  psf.FWHM_sky = fwhmSky;
  psf.FWHM_res = fwhmRes;

  
  /////////////////////
  // .... Display and verbose
  //////////////////////////////////////////////////////////
  if(disp){
    meth = "Full analytic";
    if(psfrMethod == "ls")
      meth = "TS-based LS";
    else if(psfrMethod == "mmse")
      meth = "TS-based MMSE";
    
    l = cam.nPixelsCropped * cam.pixSize;
    window,0; clr;logxy,0,0; pli, *psf.ncpa,-l/2,-l/2,l/2,l/2,cmin=0;
    pltitle,"Best bench PSF with SR = " + var2str(arrondi(100*psf.SR_ncpa,1))+"%";
    xytitles,"Arcsec","Arcsec";
  
    window,1; clr;pli, *psf.res,-l/2,-l/2,l/2,l/2,cmin=0,cmax = psf.SR_sky;
    pltitle,meth +" PSF with SR = " + var2str(arrondi(100*psf.SR_res,1))+"%";
    xytitles,"Arcsec","Arcsec";
  
    window,2; clr;pli, *psf.sky,-l/2,-l/2,l/2,l/2,cmin=0,cmax = psf.SR_sky;
    pltitle,"On-sky PSF with SR = " + var2str(arrondi(100*psf.SR_sky,1)) +"%";
    xytitles,"Arcsec","Arcsec";
  
    window,3; clr; pli, *psf.diff,-l/2,-l/2,l/2,l/2,cmin=0,cmax = psf.SR_sky;
    pltitle,"Residual of the reconstruction ";
    xytitles,"Arcsec","Arcsec";
  
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
  
    winkill,5;window,5,style="aanda.gs",dpi=90;clr;
    plg, *psf.EE_res, boxsize,color=[128,128,128];
    plg, *psf.EE_sky, boxsize;
    plg, [100,100],[-0.1,max(boxsize)*1.05],type=2,marks=0;
    xytitles,"Angular separation from center (arcsec)","Ensquared Energy (%)";
    plt,"A: "+ meth + " PSF",0.5,30,tosys=1,color=[128,128,128];
    plt,"B: On-sky PSF",0.5,25,tosys=1;
    range,0,105;
    limits,-0.1,max(boxsize)*1.05;

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

    if(strpart(rtc.aoMode,1:4) == "MOAO"){
      winkill,7;window,7,style="aanda.gs",dpi=90;clr;
      displayLayers,atm.cnh,atm.altitude,col=[char(242)],percent=1,thick=.4;
      displayLayers,-(*rtc.skyProfile)(,1),(*rtc.skyProfile)(,2),col=[char(241)],percent=1;
      limits,-100,100;
      range,-1,20;
      plt,"Calibration on-sky profile",-80,15,tosys=1;
      plt,"Post retrieved profile",30,15,tosys=1,color=[128,128,128];
    }
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
    write,format="Diff of psf         = %.3g%s\n", 100*psf.diffPsf,"%";
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
  pracResults(11) = &(100*[psf.SR_sky,psf.SR_res,budget.SRmar,budget.SRpar,budget.SRborn,psf.SR_tomo,psf.SR_fit,psf.SR_bw,psf.SR_stats,psf.SR_ncpa,psf.diffPsf]);

  if(writeRes){
    savingDir = "results/pracResults/resultsWith" + method + "_Method/";
    if(!direxist(savingDir))
      system,"mkdir " + savingDir;
    writefits, savingDir + "pracResults_" + method + "_" + strpart(procDir,1:10) + "_" + timedata + ".fits",pracResults;

  }
  
  return pracResults;
}

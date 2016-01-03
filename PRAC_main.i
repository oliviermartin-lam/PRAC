include, "PRAC_telescope.i"; 
include, "PRAC_noise.i";
include, "PRAC_dm.i"; 
include, "PRAC_bwfit.i"; 
include, "PRAC_tomography.i";
include, "PRAC_psf.i";
include, "PRAC_lmfit.i";
include, "PRAC_learn.i";
include, "PRAC_errorbreakdown.i";

include, "PRAC_tools_context.i";
include, "PRAC_tools_fits.i";
include, "PRAC_tools_math.i";
include, "PRAC_tools_turbu.i";
include, "PRAC_tools_display.i";

pldefault,font="times";
pltitle_font = "times";

/*
slopestl_2013-09-18_02h30m25s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h30m48s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h31m11s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h31m34s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h31m57s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h32m24s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h32m47s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h33m10s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h33m36s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h34m00s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h34m23s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h34m45s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h35m09s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h35m32s_script292_scao_glao4L3N_moao4l3n.fits
slopestl_2013-09-18_02h35m59s_script292_scao_glao4L3N_moao4l3n.fits

slopestl_2013-09-18_02h39m09s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h39m32s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h39m55s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h40m18s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h40m42s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h41m09s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h41m32s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h41m55s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h42m18s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h42m41s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h43m05s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h43m28s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h43m50s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h44m14s_script293_scao_moao0L3N_moao4l3n.fits
slopestl_2013-09-18_02h44m37s_script293_scao_moao0L3N_moao4l3n.fits
*/

func PRAC_main(timedata,mode=,Dir=,verb=,disp=,budgetonly=)
/*DOCUMENT p = PRAC_main("02h31m34s",verb=1,disp=0); p = PRAC_main("02h41m32s",verb=1,disp=1);
 
  to be done:
  - correlation BW and TOMO to be fixed
  - include vibration from data.learn.tracking
  - noise modelization
  - include a simulation mode
 */
{
  extern NBWFS_DEF,NL_DEF;

  tic,1;
  
  include, "PRAC_config.i",1;
  if(Dir){
    procDir = Dir;
    dataDirRoot = dataDir + procDir;
  }
  if(is_void(verb)){
    verb=0;
  }

  //...... Initializes the data struct .....//
  include, "PRAC_struct_init.i",1;
  include, "PRAC_struct_config.i",1;
  define_structs, timedata,verb=verb;
  if(data.wfs(data.its).type == 0){
    return 0;
  }

  if(!budgetonly){
    //..... Derivation of the perfect telescope OTF .....//
    OTF_tel_tmp = OTF_telescope();

    //...... Derivation of the phase structure functions DPHI_bwfit from mixed bandwidth and fitting error .....//
    //Computes an interaction matrix
    mia = intermat(data.dm);

    // Inversion of the interaction matrix -> command matrix
    mca = computeCommandMat(mia, nmf=-1, condi=30., disp=0);

    // Derivation of the phase structure function from the bw + fit errors in radians^2
    if(is_void(geo)) geo = "square";
  
    DPHI_bwfit = compute_DphiBwFitting(data.turbu.r0ir, data.turbu.L0, data.turbu.v,data.turbu.dir,data.rtc.Fe,data.rtc.delay,data.rtc.gain,data.rtc.BP,geo,mode=data.rtc.obsmode,verb=verb);

    
    //...... Derivation of the OTF from static aberration except NCPA .....//
    OTF_stats_tmp = PRAC_OTF_stats(mca,PSF_stats);
  }
  //..... Derivation of the error tomographic residual covariance matrix in arcsec^2 .....//
  Cee = computesCeeMatrix(data.rtc.obsmode,verb=verb);
  if(!budgetonly){
    //..... Derivation of the phase struture function of the tomography only in volts^2 .....//

    //computes the covariance matrix of voltages
    Cvv = propagateError2Dm(Cee, mca,verb=verb );

  
    //computing the phase structure function from the tomo error in rd^2
    DPHI_tomo = dphiFromUij(Cvv,mode,verb=verb);

    //...... Derivation of the OTF from bandwidth,fitting and tomographic error .....//
    OTF_dyn_tmp =  exp(-0.5*(DPHI_tomo + DPHI_bwfit));
   
    //manages the size of the fourier support to get the same pixel scale
    //of the IR cam as on-sky: N' = (uz*N)/uz'
    npix_sr = data.camir.npix_sr;
    npix = data.fourier.npix;
    npix2 = int(data.fourier.uz * npix/data.camir.uz);

    if(npix2 > npix){
      OTF_dyn = OTF_stats = OTF_tel = array(0.,npix2,npix2);
      nm = (npix2-npix)/2+1;
      np = (npix2+npix)/2;
      OTF_dyn(nm:np,nm:np)   = roll(OTF_dyn_tmp);
      OTF_tel(nm:np,nm:np)   = roll(OTF_tel_tmp);
      OTF_stats(nm:np,nm:np) = roll(OTF_stats_tmp);
    }else{
      OTF_dyn = OTF_stats = OTF_tel = array(0.,npix,npix);
      nm = (npix-npix2)/2+1;
      np = (npix+npix2)/2;
      OTF_dyn = roll(roll(OTF_dyn_tmp)(nm:np,nm:np));
      OTF_tel =  roll(roll(OTF_tel_tmp)(nm:np,nm:np)+1e-15);
      OTF_stats =  roll(roll(OTF_stats_tmp)(nm:np,nm:np));
    }

    OTF_ts = OTF_dyn * OTF_stats/(OTF_tel);
  
    //product of all FTO but ncpa's one
    PSF_ts = roll( fft(OTF_ts).re );
    //crops to compare to the on-sky psf
    nm = (npix2 - npix_sr)/2+1;
    np = (npix2 + npix_sr)/2;
    PSF_ts = PSF_ts(nm:np,nm:np);
    OTF_ts = fft(roll(PSF_ts)).re/numberof(PSF_ts);
  }
  //..... Gets the best PSF on bench .....//
  tmp = giveNCPADataDir(procDir);
  DirNCPA = tmp(1);
  timencpa = tmp(2);
  dataDirRoot = dataDir +  DirNCPA;
  OTF_ncpa = PRAC_OTF_ncpa(timencpa,npix_sr,SR_best,PSF_ncpa,disp=disp);

  if(!budgetonly){
    //..... Derivation of the residual OTF .....//
    OTF_res = OTF_ts * OTF_ncpa;

    //..... Derivation of the residual PSF .....//
    PSF_res  = roll( fft(OTF_res).re );
    PSF_res /= sum(PSF_res);
    airyPeak = (pi/4)*(1 - data.tel.obs*data.tel.obs) * data.camir.uld * data.camir.uld ;
    PSF_res /= airyPeak;
  }
    //diff
    dataDirRoot = dataDir +  procDir;
    suff_ir = readFitsKey(pathdata,"IRFILE");
    restorefits,"ir",suff_ir,pathir,fake=1;
    PSF_sky = processOnePSF(suff_ir,SR_sky,snrsky,box = npix_sr,disp=disp);
    data.uncertainties.srsky = snrsky;
  

  //.... Error breakdown .....//
  data.budget.ncpa = sqrt(- (data.camir.lambda_ir*1e9/2/pi)^2 * log(SR_best));
  if(verb) write,"Computing the error breakdown";
  PRAC_errorbreakdown,verb=verb;

  data.budget.SRsky = 100*SR_sky;
  if(!budgetonly){
    SR_tomo = exp(-(2*pi*dphi2rms(DPHI_tomo) *1e-9/data.camir.lambda_ir)^2);
    SR_bwfit = exp(-(2*pi*dphi2rms(DPHI_bwfit) *1e-9/data.camir.lambda_ir)^2);
    SR_stats = max(PSF_stats/sum(PSF_stats)/airyPeak);
    data.budget.SRres = SR_res = max(PSF_res);
    PSF_diff = abs(PSF_sky - PSF_res);
    diffPSF = sum(PSF_diff(*)^2)/sum(PSF_sky(*)^2);
  
 
    if(disp){
   
      l = npix_sr * data.camir.uz;
      window,0; clr; palette,"stern.gp";logxy,0,0; pli, PSF_ncpa,-l/2,-l/2,l/2,l/2,cmin=0;
      pltitle,"Best bench PSF with SR = " + var2str(arrondi(100*SR_best,1))+"%";
      xytitles,"Arcsec","Arcsec";
  
      window,1; clr;palette,"stern.gp";pli, PSF_res,-l/2,-l/2,l/2,l/2,cmin=0,cmax = SR_sky;
      pltitle,"Reconstructed PSF with SR = " + var2str(arrondi(100*SR_res,1))+"%";
      xytitles,"Arcsec","Arcsec";
    
      window,2; clr;palette,"stern.gp"; pli, PSF_sky,-l/2,-l/2,l/2,l/2,cmin=0,cmax = SR_sky;
      pltitle,"On-sky PSF with SR = " + var2str(arrondi(100*SR_sky,1)) +"% +/- "+ var2str(arrondi(100*snrsky,1)) + "%";
      xytitles,"Arcsec","Arcsec";

      window,3; clr;palette,"stern.gp"; pli, PSF_diff,-l/2,-l/2,l/2,l/2,cmin=0,cmax = SR_sky;
      pltitle,"Residual of the reconstruction ";
      xytitles,"Arcsec","Arcsec";

      winkill,4;window,4,style="aanda.gs",dpi=90;clr;
      y = [data.budget.res,data.budget.tomo,data.budget.alias,data.budget.noise,data.budget.bw,data.budget.fit,data.budget.ol,data.budget.static,data.budget.ncpa];
      labs = ["!s_IR","!s_tomo","!s_alias","!s_noise","!s_bw","!s_fit","!s_ol","!s_static","!s_ncpa"];
      plotsBarDiagram,y,labs,col1=[char(241)],title=1;
    }

    boxsize = EE = EE_sky = span(1,npix_sr-3,npix_sr) * data.camir.uz;
    for(i=1; i<=npix_sr; i++){ EE(i) = getEE( 100*PSF_res/sum(PSF_res), data.camir.uz, boxsize(i) );}
    for(i=1; i<=npix_sr; i++){ EE_sky(i) = getEE(100*PSF_sky/sum(PSF_sky), data.camir.uz, boxsize(i) );}
    if(disp){
      winkill,5;window,5,style="aanda.gs",dpi=90;clr; palette,"earth.gp";
      plg, EE, boxsize,color=[128,128,128];
      plg, EE_sky, boxsize;
      xytitles,"Size of the box (arcsec)","Ensquared Energy (%)";
      plt,"A: Reconstructed PSF",0.5,30,tosys=1,color=[128,128,128];
      plt,"B: On-sky PSF",0.5,25,tosys=1;
      range,0,105;
      limits,-0.1,max(boxsize)*1.05;
    }
    tf = tac(1);
    if(verb){
      write,format="SR on-sky           = %.3g%s +/- %.2g%s\n", 100*SR_sky,"%",100*snrsky,"%";
      write,format="SR reconstructed    = %.3g%s\n", 100*SR_res,"%";
      write,format="SR Mar. from budget = %.3g%s\n", data.budget.SRmar,"%";
      write,format="SR Par. from budget = %.3g%s\n", data.budget.SRpar,"%";
      write,format="SR Bor. from budget = %.3g%s\n", data.budget.SRborn,"%";
      write,format="SR bwfit            = %.3g%s\n", 100*SR_bwfit,"%";
      write,format="SR tomo+alias+noise = %.3g%s\n", 100*SR_tomo,"%";
      write,format="SR static           = %.3g%s\n", 100*SR_stats,"%";
      write,format="SR ncpa             = %.3g%s\n", 100*SR_best,"%";
      write,format="Diff of psf         = %.3g%s\n", 100*diffPSF,"%";
      write,format="Residual error      = %.4g nm rms\n", data.budget.res;
      write,format="Tomographic error   = %.4g nm rms\n", data.budget.tomo;
      write,format="Aliasing error      = %.4g nm rms\n", data.budget.alias;
      write,format="Noise error         = %.4g nm rms\n", data.budget.noise;
      write,format="Bandwidth error     = %.4g nm rms\n", data.budget.bw;
      write,format="Fitting error       = %.4g nm rms\n", data.budget.fit;
      write,format="Go-to error         = %.4g nm rms\n", data.budget.ol;
      write,format="Static error        = %.4g nm rms\n", data.budget.static;
      write,format="NCPA error          = %.4g nm rms\n", data.budget.ncpa;
      write,format="PSF reconstruction done on %.3g s\n",tf;
    }

    //merging results
    PRAC_res = array(pointer,7);
    //Data identity
    PRAC_res(1) = &[timedata,data.rtc.aomode,data.rtc.obsmode,data.rtc.rectype];
    //global parameters
    PRAC_res(2) = &[data.turbu.r0vis,data.turbu.L0,data.turbu.v,tf];
    //turbulence profile
    PRAC_res(3) = &[data.learn.cnh,data.learn.altitude,data.learn.l0h];
    //tracking
    PRAC_res(4) = &data.learn.tracking;
    //error budget
    PRAC_res(5) = &[data.budget.res,data.budget.fit,data.budget.bw,data.budget.tomo,data.budget.noise,data.budget.static,data.budget.ncpa];
    //Strehl ratios
    PRAC_res(6) = &(100*[SR_sky,snrsky,SR_res,data.budget.SRmar,data.budget.SRpar,data.budget.SRborn,SR_tomo,SR_bwfit,SR_stats,SR_best,diffPSF]);
    //Ensquared Energy
    PRAC_res(7) = &avg(abs(100*(EE-EE_sky)/EE_sky));
    
    return PRAC_res;
  }else{
    p = array(pointer,6);
    p(1) = &data.rtc.aomode;
    p(2) = &[data.turbu.r0vis,data.turbu.L0,data.turbu.v];
    p(3) = &[data.learn.cnh,data.learn.altitude,data.learn.l0h,data.learn.vh];
    p(4) = &[data.budget.res,data.budget.fit,data.budget.bw,data.budget.tomo,data.budget.alias,data.budget.noise,data.budget.ol,data.budget.static,data.budget.ncpa,data.budget.SRsky,data.budget.SRmar,data.budget.SRpar,data.budget.SRborn,data.budget.vib];
    p(5) = &PSF_sky;
    p(6) = &data.budget.ved;
    return p;
  }
}





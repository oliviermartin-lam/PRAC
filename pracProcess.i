include, "pracMain.i";
                  

func comparePSFRmethods(timedata,Dir=,disp=)
/* DOCUMENT comparePSFRmethods,timedata,Dir="2013_09_13_onsky/"

 */
{
  res_ts  = pracMain(timedata,Dir=Dir,psfrMethod="ts",disp=disp);
  res_est = pracMain(timedata,Dir=Dir,psfrMethod="estimation",disp=disp);
  res_rtc = pracMain(timedata,Dir=Dir,psfrMethod="rtc",disp=disp);

  ///////////
  // Relative Error on EE
  ////////////////////////////////////////////////
  
  EEsky  = (*res_ts(7))(,1);
  EEts   = (*res_ts(7))(,2);
  EEest  = (*res_est(7))(,2);
  EErtc  = (*res_rtc(7))(,2);
  EEtsdiff = 100*abs(EEts-EEsky)/EEsky;
  EEestdiff = 100*abs(EEest-EEsky)/EEsky;
  EErtcdiff = 100*abs(EErtc-EEsky)/EEsky;
  
  winkill,10;window,10,dpi=90,style="aanda.gs";
  logxy,1,0;gridxy,1,1;
  boxwidth = 1e3*span(1,cam.nPixelsCropped-3.,cam.nPixelsCropped) * cam.pixSize;
  plg,EEtsdiff ,boxwidth;
  plg,EEestdiff,boxwidth,color = [128,128,128];
  plg,EErtcdiff,boxwidth,color = [100,100,100];
  xytitles,"Integrated box width [mas]","Relative error on EE [%]";
  m = max(max(EEtsdiff),max(EEestdiff),max(EErtcdiff));
  plt,"A: TS-based  PSF",.5e3,.9*m,tosys=1;
  plt,"B: Estimated PSF",.5e3,.8*m,tosys=1,color = [128,128,128];
  plt,"C: RTC-based PSF",.5e3,.7*m,tosys=1,color = [100,100,100];
  pltitle,rtc.aoMode + " at " + timedata;
  pdf,"results/EEerror_" + timedata;

  ///////////
  //    OTF
  ////////////////////////////////////////////////
  
  os   = (*res_ts(8))(,,1);
  ots  = (*res_ts(8))(,,2);
  oest = (*res_est(8))(,,2);
  ortc = (*res_rtc(8))(,,2);
   
  os_avg   = circularAveragePsf(os);os_avg/=max(os_avg);
  ots_avg  = circularAveragePsf(ots);ots_avg/=max(ots_avg);
  oest_avg = circularAveragePsf(oest);oest_avg/=max(oest_avg);
  ortc_avg = circularAveragePsf(ortc);ortc_avg/=max(ortc_avg);

  winkill,11;window,11,style="aanda.gs",dpi=90;clr;
  dl =  (indgen(cam.nPixelsCropped/2) * tel.foV/tel.nPixels);
  n = numberof(dl);gridxy,1,1;
  plg, 10*log10(abs(os_avg)),dl,type=2;
  plg, 10*log10(abs(ots_avg)),dl;
  plg, 10*log10(abs(oest_avg)),dl;
  plg, 10*log10(abs(ortc_avg)),dl;
  xytitles,"Normalized frequency [D/!l units]","Circularly averaged OTF [dB]";
  limits,-.1,1.;range,-30,.1;
  plt,"A: Sky OTF",0,-10,tosys=1;
  plt,"B: TS-based  OTF",0,-14,tosys=1;
  plt,"C: Estimated OTF",0,-16,tosys=1,color = [128,128,128];
  plt,"D: RTC-based OTF",0,-18,tosys=1,color = [100,100,100];
  pltitle,rtc.aoMode + " at " + timedata;
  
  pdf,"results/OTFerror_" + timedata;

  ///////////
  //    PSF
  ////////////////////////////////////////////////
  
  os   = (*res_ts(6))(,,1);
  ots  = (*res_ts(6))(,,2);
  oest = (*res_est(6))(,,2);
  ortc = (*res_rtc(6))(,,2);
   
  //os_avg   = circularAveragePsf(os);//os_avg/=max(os_avg);
  ots_avg  = circularAveragePsf(ots-os);//ots_avg/=max(ots_avg);
  oest_avg = circularAveragePsf(oest-os);//oest_avg/=max(oest_avg);
  ortc_avg = circularAveragePsf(ortc-os);//ortc_avg/=max(ortc_avg);

  winkill,13;window,13,style="aanda.gs",dpi=90;clr;
  dl     = span(1,cam.nPixelsCropped/2,cam.nPixelsCropped/2) * cam.pixSize;
  n      = numberof(dl);gridxy,1,1;
  difts  = 10*log10(abs(ots_avg));
  difest = 10*log10(abs(oest_avg));
  difrtc = 10*log10(abs(ortc_avg));
  plg, difts, dl,type=2;
  plg, difest,dl,color = [128,128,128];
  plg, difrtc,dl,color = [100,100,100];
  xytitles,"Angular distance from center [arcsec]","Reconstruction absolute residue [dB]";
  w = where(dl<=.8);
  m = max(max(difts(w)),max(difest(w)),max(difrtc(w)));
  limits,-.1,.8;//range,0,1.05*m;
  //plt,"A: Sky PSF (dashed line)",0.4,-5,tosys=1;
  plt,"A: TS-based  PSF (dashed line)",0.4,-25,tosys=1;
  plt,"B: Estimated PSF",0.4,-30,tosys=1,color = [128,128,128];
  plt,"C: RTC-based PSF",0.4,-35,tosys=1,color = [100,100,100];
  pltitle,rtc.aoMode + " at " + timedata;
  
  pdf,"results/PSFerror_" + timedata;
  

  ///////////
  // Scalar results on performance
  ////////////////////////////////////////////////

  write,format="r0              = %.3g cm\n", 100*(*res_ts(2))(1);
  write,format="L0              = %.3g m\n", (*res_ts(2))(2);
  write,format="v               = %.3g m/s\n", (*res_ts(2))(3);
  
  write,format="Sky  Strehl     = %.3g%s\n", 100.*(*res_ts(11))(1)," %";
  write,format="TS   Strehl     = %.3g%s\n", 100.*(*res_ts(11))(2)," %";
  write,format="Est. Strehl     = %.3g%s\n", 100.*(*res_est(11))(2)," %";
  write,format="RTC  Strehl     = %.3g%s\n", 100.*(*res_rtc(11))(2)," %";

  write,format="Sky  FWHM       = %.4g mas\n", 1e3*(*res_ts(11))(11);
  write,format="TS   FWHM       = %.4g mas\n", 1e3*(*res_ts(11))(12);
  write,format="Est. FWHM       = %.4g mas\n", 1e3*(*res_est(11))(12);
  write,format="RTC  FWHM       = %.4g mas\n", 1e3*(*res_rtc(11))(12);

  write,format="TS   Chi2       = %.4g\n", (*res_ts(11))(13);
  write,format="Est. Chi2       = %.4g\n", (*res_est(11))(13);
  write,format="RTC  Chi2       = %.4g\n", (*res_rtc(11))(13);

  write,format="Sky  Mof. prof. = [%.3g, %.3g, %.3g,%.3g,%.3g]\n",(*res_ts(12))(1,1),(*res_ts(12))(2,1),(*res_ts(12))(3,1),(*res_ts(12))(4,1),(*res_ts(12))(5,1);
  write,format="TS   Mof. prof. = [%.3g, %.3g, %.3g,%.3g,%.3g]\n", (*res_ts(12))(1,2),(*res_ts(12))(2,2),(*res_ts(12))(3,2),(*res_ts(12))(4,2),(*res_ts(12))(5,2);
  write,format="Est. Mof. prof. = [%.3g, %.3g, %.3g,%.3g,%.3g]\n",(*res_est(12))(1,2),(*res_est(12))(2,2),(*res_est(12))(3,2),(*res_est(12))(4,2),(*res_est(12))(5,2);
  write,format="RTC  Mof. prof. = [%.3g, %.3g, %.3g,%.3g,%.3g]\n", (*res_rtc(12))(1,2),(*res_rtc(12))(2,2),(*res_rtc(12))(3,2),(*res_rtc(12))(4,2),(*res_rtc(12))(5,2);

  ////////
  // Moffat profile
  ////////////////////////////////////////////////

  winkill,12;window,12,style="aanda.gs",dpi=90;clr;
  gridxy,1,1;
  y1  = (*res_ts(12))(1:-1,1);
  y2  = (*res_ts(12))(1:-1,2);
  y3  = (*res_est(12))(1:-1,2);
  y4  = (*res_rtc(12))(1:-1,2);
  tab = [100*abs((y2-y1)/y1),100*abs((y3-y1)/y1),100*abs((y4-y1)/y1)];
  multiBarDiagram,tab,["I_0","!a_x","!a_y","!b"];
  xytitles," ","Relative error on Moffat parameters [%]";
  //plt,"Sky PSF",.2,max(tab),tosys=1,color=[192,192,192];
  plt,"TS-based PSF",.2,max(tab)*.9,tosys=1,color=[128,128,128];
  plt,"Estimated PSF",.2,max(tab)*.8,tosys=1,color=[96,96,96];
  plt,"RTC-based PSF",.2,max(tab)*.7,tosys=1,color=[64,64,64];
  pltitle,rtc.aoMode + " at " + timedata;
  pdf,"results/Moffat_" + timedata;

}
 
/*
 __  __    _    ____ ____ _____     _______ 
|  \/  |  / \  / ___/ ___|_ _\ \   / / ____|
| |\/| | / _ \ \___ \___ \| | \ \ / /|  _|  
| |  | |/ ___ \ ___) |__) | |  \ V / | |___ 
|_|  |_/_/   \_\____/____/___|  \_/  |_____|
                                            
 _        _   _   _ _   _  ____ _   _ ___ _   _  ____ 
| |      / \ | | | | \ | |/ ___| | | |_ _| \ | |/ ___|
| |     / _ \| | | |  \| | |   | |_| || ||  \| | |  _ 
| |___ / ___ \ |_| | |\  | |___|  _  || || |\  | |_| |
|_____/_/   \_\___/|_| \_|\____|_| |_|___|_| \_|\____|
                                                      
*/

func pracAll(method)
/* DOCUMENT pracAll,"rtc";

 */
{
  
  include,"pracConfig.i",1;
  //Dir  = listFile(dataDir);
  //Dir  = Dir(sort(Dir));
  //  Dir  = Dir(17:21);
  Dir = ["2013_09_13_onsky","2013_09_15_onsky","2013_09_16_onsky","2013_09_17_onsky","2013_09_18_onsky"];

  ndir = numberof(Dir);
  
  goodDir = "results/pracResults/";
  if(!direxist(goodDir)) system,"mkdir " + goodDir;
  
  for(i=1;i<=ndir;i++){
    procDir = Dir(i);
    dataDirRoot = dataDir + procDir;
    //takes all slopestl files
    pathstl = listVersion(dataDirRoot,"fits","slopestl");
    if(is_array(pathstl)){
        //keeps only the script files
        w = where(strfind("script",pathstl)(0,) != -1);
        pathstl = pathstl(w);
        ntl = numberof(pathstl);
        for(j=1;j<=ntl;j++){
          //extracts the date
          timetl = extractTime(pathstl(j));
          mode   = strcase(1,readFitsKey(pathstl(j),"OBS_MODE"));
          p = pracMain(timetl,Dir=procDir,verb=0,disp=0,psfrMethod=method,writeRes=1);
          write,format="%s","Job done: " + var2str(100.*j/ntl) + "%\r";
        }
    }
  }
}

/*
 ____ _____  _  _____ ___ ____ _____ ___ ____ ____  
/ ___|_   _|/ \|_   _|_ _/ ___|_   _|_ _/ ___/ ___| 
\___ \ | | / _ \ | |  | |\___ \ | |  | | |   \___ \ 
 ___) || |/ ___ \| |  | | ___) || |  | | |___ ___) |
|____/ |_/_/   \_\_| |___|____/ |_| |___\____|____/ 
                                                    
*/

func psfrHistogram(void)
/* DOCUMENT
 */
{
  include,"pracConfig.i",1;
  savingDir  = "results/pracResults/resultsWithestimation_Vii_Method/";
  list       = listFile(savingDir);
  nfile      = numberof(list);  
  obsmode = SR_sky = SR_est = [];
  chiest = chits = chirtc = [];
  I0est = axest = ayest = best = thest =  [];
  I0sky = axsky = aysky = bsky = thsky =  [];
  
  for(i=1;i<=nfile;i++){
    //loading results
    p = readfits(savingDir + list(i));
    //grabbing aomode
    mode = strchar(*p(1))(2);
    obsmode = grow(obsmode,mode);
    //grabbing the Strehl ratios
    SR_sky = grow(SR_sky,(*p(11))(1));
    SR_est = grow(SR_est,(*p(11))(2)) ;
    chiest = grow(chiest,(*p(11))(13)) ;
    I0sky    = grow(I0sky,(*p(12))(1,1));
    I0est    = grow(I0est,(*p(12))(1,2));
    axsky     = grow(axsky,(*p(12))(2,1));
    axest     = grow(axest,(*p(12))(2,2));
    aysky     = grow(aysky,(*p(12))(3,1));
    ayest     = grow(ayest,(*p(12))(3,2));
    bsky      = grow(bsky,(*p(12))(4,1));
    best      = grow(best,(*p(12))(4,2));
    thsky     = grow(thsky,(*p(12))(5,1));
    thest     = grow(thest,(*p(12))(5,2));
    write,format="Files grabbed :%.3g%s\r",100.*i/nfile,"%";
  }

  savingDir  = "results/pracResults/resultsWithts_Vii_Method/";
  list       = listFile(savingDir);
  nfile      = numberof(list);  
  SR_ts= I0ts = axts =  ayts = bts = thts =  [];
  for(i=1;i<=nfile;i++){
    //loading results
    p = readfits(savingDir + list(i));
    SR_ts = grow(SR_ts,(*p(11))(2)) ;
    chits = grow(chits,(*p(11))(13)) ;
    I0ts    = grow(I0ts,(*p(12))(1,2));
    axts     = grow(axts,(*p(12))(2,2));
    ayts     = grow(ayts,(*p(12))(3,2));
    bts      = grow(bts,(*p(12))(4,2));
    thts     = grow(thts,(*p(12))(5,2));
    write,format="Files grabbed :%.3g%s\r",100.*i/nfile,"%";
  }

  savingDir  = "results/pracResults/resultsWithrtc_Vii_Method/";
  list       = listFile(savingDir);
  nfile      = numberof(list);  
  SR_rtc = I0rtc = axrtc = ayrtc = brtc = thrtc =  [];
  for(i=1;i<=nfile;i++){
    //loading results
    p = readfits(savingDir + list(i));
    SR_rtc = grow(SR_rtc,(*p(11))(2)) ;
    chirtc = grow(chirtc,(*p(11))(13)) ;
    I0rtc    = grow(I0rtc,(*p(12))(1,2));
    axrtc     = grow(axrtc,(*p(12))(2,2));
    ayrtc     = grow(ayrtc,(*p(12))(3,2));
    brtc      = grow(brtc,(*p(12))(4,2));
    thrtc     = grow(thrtc,(*p(12))(5,2));
    write,format="Files grabbed :%.3g%s\r",100.*i/nfile,"%";
  }

  ///////
  // Strehl
  ////////////////////////////////////////
  
  nc = 5;
  difSRts  = 100*(SR_ts -SR_sky)/SR_sky;
  w = where(chits <nc & abs(difSRts) < 100);
  difSRts =  difSRts(w);
  difSRest = 100*(SR_est-SR_sky)/SR_sky;
  w = where(chiest <nc & abs(difSRest) < 100);
  difSRest =  difSRest(w);
  difSRrtc = 100*(SR_rtc-SR_sky)/SR_sky;
  w = where(chirtc <nc & abs(difSRrtc) < 100);
  difSRrtc =  difSRrtc(w);
  // Histograms
  n = 20;
  winkill,1;window,1,dpi=90,style="aanda.gs";gridxy,1,1;
  histo,difSRts,n;
  histo,difSRest,n,color="blue";
  histo,difSRrtc,n,color="red";
  xytitles,"Total error on Strehl [%]","Counts";
  m = limits();
  plt,"A: TS-based method",.45*m(2),.9*m(4),tosys=1;
  plt,"B: Estimation method",.45*m(2),.8*m(4),color="blue",tosys=1;
  plt,"C: RTC-based method",.45*m(2),.7*m(4),color="red",tosys=1;
  pdf,"results/statistics/strehlHistogram";
  
  ///////
  // Moffat I0
  ////////////////////////////////////////
  
  difI0ts  = 100*(I0ts -I0sky)/I0sky;
  w = where(chits <nc & abs(difI0ts) < 100);
  difI0ts =  difI0ts(w);
  difI0est = 100*(I0est-I0sky)/I0sky;
  w = where(chiest <nc & abs(difI0est) < 100);
  difI0est =  difI0est(w);
  difI0rtc = 100*(I0rtc-I0sky)/I0sky;
  w = where(chirtc <nc & abs(difI0rtc) < 100);
  difI0rtc =  difI0rtc(w);
  // Histograms
  n = 20;
  winkill,2;window,2,dpi=90,style="aanda.gs";gridxy,1,1;
  histo,difI0ts,n;
  histo,difI0est,n,color="blue";
  histo,difI0rtc,n,color="red";
  xytitles,"Total error on I_0_ [%]","Counts";
  m = limits();
  plt,"A: TS-based method",.45*m(2),.9*m(4),tosys=1;
  plt,"B: Estimation method",.45*m(2),.8*m(4),color="blue",tosys=1;
  plt,"C: RTC-based method",.45*m(2),.7*m(4),color="red",tosys=1;
  pdf,"results/statistics/I0Histogram";

  ///////
  // Moffat alpha x
  ////////////////////////////////////////

  difaxts  = 100*(axts -axsky)/axsky;
  w = where(chits <nc & abs(difaxts) < 100);
  difaxts =  difaxts(w);
  difaxest = 100*(axest-axsky)/axsky;
  w = where(chiest <nc & abs(difaxest) < 100);
  difaxest =  difaxest(w);
  difaxrtc = 100*(axrtc-axsky)/axsky;
  w = where(chirtc <nc & abs(difaxrtc) < 100);
  difaxrtc =  difaxrtc(w);
  // Histograms
  n = 20;
  winkill,3;window,3,dpi=90,style="aanda.gs";gridxy,1,1;
  histo,difaxts,n;
  histo,difaxest,n,color="blue";
  histo,difaxrtc,n,color="red";
  xytitles,"Total error on !a_x [%]","Counts";
  m = limits();
  plt,"A: TS-based method",.45*m(2),.9*m(4),tosys=1;
  plt,"B: Estimation method",.45*m(2),.8*m(4),color="blue",tosys=1;
  plt,"C: RTC-based method",.45*m(2),.7*m(4),color="red",tosys=1;
  pdf,"results/statistics/alphaHistogram";

  ///////
  // Moffat alpha y
  ////////////////////////////////////////

  difayts  = 100*(ayts -aysky)/aysky;
  w = where(chits <nc & abs(difayts) < 100);
  difayts =  difayts(w);
  difayest = 100*(ayest-aysky)/aysky;
  w = where(chiest <nc & abs(difayest) < 100);
  difayest =  difayest(w);
  difayrtc = 100*(ayrtc-aysky)/aysky;
  w = where(chirtc <nc & abs(difayrtc) < 100);
  difayrtc =  difayrtc(w);
  // Histograms
  n = 20;
  winkill,4;window,4,dpi=90,style="aanda.gs";gridxy,1,1;
  histo,difayts,n;
  histo,difayest,n,color="blue";
  histo,difayrtc,n,color="red";
  xytitles,"Total error on !a_y [%]","Counts";
  m = limits();
  plt,"A: TS-based method",.45*m(2),.9*m(4),tosys=1;
  plt,"B: Estimation method",.45*m(2),.8*m(4),color="blue",tosys=1;
  plt,"C: RTC-based method",.45*m(2),.7*m(4),color="red",tosys=1;
  pdf,"results/statistics/alphaYHistogram";
  

  ///////
  // Moffat beta
  ////////////////////////////////////////

  difbts  = 100*(bts -bsky)/bsky;
  w = where(chits <nc & abs(difbts) < 100);
  difbts =  difbts(w);
  difbest = 100*(best-bsky)/bsky;
  w = where(chiest <nc & abs(difbest) < 100);
  difbest =  difbest(w);
  difbrtc = 100*(brtc-bsky)/bsky;
  w = where(chirtc <nc & abs(difbrtc) < 100);
  difbrtc =  difbrtc(w);
  // Histograms
  n = 20;
  winkill,5;window,5,dpi=90,style="aanda.gs";gridxy,1,1;
  histo,difbts,n;
  histo,difbest,n,color="blue";
  histo,difbrtc,n,color="red";
  xytitles,"Total error on !b [%]","Counts";
  m = limits();
  plt,"A: TS-based method",.45*m(2),.9*m(4),tosys=1;
  plt,"B: Estimation method",.45*m(2),.8*m(4),color="blue",tosys=1;
  plt,"C: RTC-based method",.45*m(2),.7*m(4),color="red",tosys=1;
  pdf,"results/statistics/betaHistogram";

  /*
  ///////
  // Moffat theta
  ////////////////////////////////////////

  difthts  = 100*(thts -thsky)/thsky;
  w = where(chits <nc & abs(difthts) < 200);
  difthts =  difthts(w);
  difthest = 100*(thest-thsky)/thsky;
  w = where(chiest <nc & abs(difthest) < 200);
  difthest =  difthest(w);
  difthrtc = 100*(thrtc-thsky)/thsky;
  w = where(chirtc <nc & abs(difthrtc) < 200);
  difthrtc =  difthrtc(w);
  // Histograms
  n = 15;
  winkill,6;window,6,dpi=90,style="aanda.gs";
  histo,difthts,n;
  histo,difthest,n,color="blue";
  histo,difthrtc,n,color="red";
  xytitles,"Total error on !q [%]","Counts";
  pdf,"results/statistics/thetaHistogram";
  */

  ///////////////////
  // STATISTICS ON RECONSTRUCTION
  ////////////////////////////////////////////////

  write,format="Average on RTC-based SR diff: %.3g\n", difSRrtc(avg);
  write,format="1-sigma on RTC-based SR diff: %.3g\n", difSRrtc(rms);
  write,format="Correlation on RTC-based SR : %.3g\n", getCorrelationFactor(SR_rtc,SR_sky)
  write,format="Average on Estimated SR diff: %.3g\n", difSRest(avg);
  write,format="1-sigma on Estimated SR diff: %.3g\n", difSRest(rms);
  write,format="Correlation on Estimated SR : %.3g\n", getCorrelationFactor(SR_ts,SR_sky)
  write,format="Average on TS-based  SR diff: %.3g\n", difSRts(avg); 
  write,format="1-sigma on TS-based  SR diff: %.3g\n", difSRts(rms);
  write,format="Correlation on TS-based SR  : %.3g\n \n", getCorrelationFactor(SR_est,SR_sky)
                    

  
  write,format="Average on RTC-based I0 diff: %.3g\n", difI0rtc(avg);
  write,format="1-sigma on RTC-based I0 diff: %.3g\n", difI0rtc(rms);
  write,format="Correlation on RTC-based I0 : %.3g\n", getCorrelationFactor(I0rtc,I0sky)
  write,format="Average on Estimated I0 diff: %.3g\n", difI0est(avg);
  write,format="1-sigma on Estimated I0 diff: %.3g\n", difI0est(rms);
  write,format="Correlation on Estimated I0 : %.3g\n", getCorrelationFactor(I0est,I0sky)
  write,format="Average on TS-based  I0 diff: %.3g\n", difI0ts(avg);  
  write,format="1-sigma on TS-based  I0 diff: %.3g\n", difI0ts(rms);
  write,format="Correlation on TS-based I0 : %.3g\n \n", getCorrelationFactor(I0ts,I0sky)

  write,format="Average on RTC-based ax diff: %.3g\n", difaxrtc(avg);
  write,format="1-sigma on RTC-based ax diff: %.3g\n", difaxrtc(rms);
  write,format="Correlation on RTC-based ax : %.3g\n", getCorrelationFactor(axrtc,axsky)
  write,format="Average on Estimated ax diff: %.3g\n", difaxest(avg);
  write,format="1-sigma on Estimated ax diff: %.3g\n", difaxest(rms);
  write,format="Correlation on Estimated ax : %.3g\n", getCorrelationFactor(axest,axsky)
  write,format="Average on TS-based  ax diff: %.3g\n", difaxts(avg);
  write,format="1-sigma on TS-based  ax diff: %.3g\n", difaxts(rms);
  write,format="Correlation on TS-based ax : %.3g\n \n", getCorrelationFactor(axts,axsky)

  write,format="Average on RTC-based ay diff: %.3g\n", difayrtc(avg);
  write,format="1-sigma on RTC-based ay diff: %.3g\n", difayrtc(rms);
  write,format="Correlation on RTC-based ay : %.3g\n", getCorrelationFactor(ayrtc,aysky)
  write,format="Average on Estimated ay diff: %.3g\n", difayest(avg);
  write,format="1-sigma on Estimated ay diff: %.3g\n", difayest(rms);
  write,format="Correlation on Estimated ay : %.3g\n", getCorrelationFactor(ayest,aysky)
  write,format="Average on TS-based  ay diff: %.3g\n", difayts(avg);
  write,format="1-sigma on TS-based  ay diff: %.3g\n", difayts(rms);
  write,format="Correlation on TS-based ay : %.3g\n \n", getCorrelationFactor(ayts,aysky)


  write,format="Average on RTC-based be diff: %.3g\n", difbrtc(avg);
  write,format="1-sigma on RTC-based be diff: %.3g\n", difbrtc(rms);
  write,format="Correlation on RTC-based be : %.3g\n", getCorrelationFactor(brtc,bsky)
  write,format="Average on Estimated be diff: %.3g\n", difbest(avg);
  write,format="1-sigma on Estimated be diff: %.3g\n", difbest(rms);
  write,format="Correlation on Estimated be : %.3g\n", getCorrelationFactor(best,bsky)
  write,format="Average on TS-based  be diff: %.3g\n", difbts(avg);
  write,format="1-sigma on TS-based  be diff: %.3g\n", difbts(rms);
  write,format="Correlation on TS-based be : %.3g\n", getCorrelationFactor(bts,bsky)

  
}

func statisticsOnPSFR(method,aomode,all=)
/* DOCUMENT statisticsOnPSFR,"ts",[],all=1
 */
{

  include,"pracConfig.i",1;

  savingDir  = "results/pracResults/resultsWith" + method + "_Vii_Method/"
  list       = listFile(savingDir);
  nfile      = numberof(list);
  pixSize    = 0.0297949;
  //pixSize   *= 128./512.;
  
  obsmode = FWHM_sky = FWHM_res = chi2 = r0 =[];
  I0sky = I0res = ares = asky = bres=bsky = [];
  SR_sky = SR_res = [];
  os = or = array(0.,512,512);
  N  = 128;
  xi = 256-N/2+1;
  xf = 256 + N/2;
    
  for(i=1;i<=nfile;i++){
    //loading results
    p = readfits(savingDir + list(i));
    //grabbing aomode
    mode = strpart(strchar(*p(1))(2),1:4);

    if(mode == aomode || all == 1){
      obsmode = grow(obsmode,mode);
      SR_sky = grow(SR_sky,(*p(11))(1));
      SR_res = grow(SR_res,(*p(11))(2)) ;
      FWHM_sky = grow(FWHM_sky,(*p(11))(11));
      FWHM_res = grow(FWHM_res,(*p(11))(12)) ;
      chi2     = grow(chi2,(*p(11))(13)) ;
      r0       = grow(r0,(*p(2))(1));
      I0sky    = grow(I0sky,(*p(12))(1,1));
      I0res    = grow(I0res,(*p(12))(1,2));
      asky     = grow(asky,(*p(12))(2,1));
      ares     = grow(ares,(*p(12))(2,2));
      bsky     = grow(bsky,(*p(12))(3,1));
      bres     = grow(bres,(*p(12))(3,2));
                      
      write,format="Files grabbed :%.3g%s\r",100.*i/nfile,"%";
    }
  }
  nc       = 50;
  w        = where(chi2 <=nc);
  obsmode  = obsmode(w);
  SR_sky   = 100*SR_sky(w);
  SR_res   = 100*SR_res(w);
  FWHM_sky = FWHM_sky(w);
  FWHM_res = FWHM_res(w);
  chi2     = chi2(w);  
  r0       = r0(w);
  I0sky    = I0sky(w);
  I0res    = I0res(w);
  asky     = asky(w);
  ares     = ares(w);
  bsky     = bsky(w);
  bres     = bres(w);

  a = sort(SR_sky);
  SR_sky = SR_sky(a);
  SR_res = SR_res(a);
  //mode = sortLabel(obsmode);
  mode = ["SCAO","MOAO","GLAO"];
  c  = ["black","red","blue"];
  
  winkill,0;window,0,dpi=90,style="aanda.gs";
  m = max(max(SR_res),max(SR_sky));
 
  for(j=1;j<=numberof(mode);j++){
    w = where(strpart(obsmode,1:4) == mode(j));
    ns       = numberof(w);
    if(is_array(w)){
      coeff = proportionalRegress(SR_res(w),SR_sky(w));
      rho   = getCorrelationFactor(SR_res(w),SR_sky(w));
      plmk, SR_res(w),SR_sky(w),marker = j,msize=.2,width = 10,color=c(j);
     
      plt,mode(j) + " ("+var2str(ns) +" samples): " +var2str(arrondi(100*rho,1)) + "% of correlation",2,(.9-j/10.)*m,color=c(j),tosys=1;
      plg,coeff*SR_sky,SR_sky,marks=0,type=2,width=4,color=c(j);
    }
  }
  
  xytitles,"H-band Sky Strehl ratio [%]","Reconstructed Strehl ratio [%]";
  gridxy,1,1;
  limits,0,m*1.05;
  range,0,m*1.05;

  meth = "Estimation method";
  if(method == "ts")
    meth = "TS-based method";
  else if(method == "rtc")
    meth = "RTC-based method";
    
  if(is_string(aomode)){
    pltitle, meth + " in " +  aomode;
    pdf,"results/srReconstructionIn" + aomode + "_" + method;
  }else{
    pltitle, meth;
    pdf,"results/statistics/srReconstructionInAllAoModes" + "_" + method;
  }
  
  
  //////
  // FWHM
  ///////////////////////////
  
  FWHM_sky = FWHM_sky*1e3;
  FWHM_res = FWHM_res*1e3;
  
  w = where(FWHM_sky < 200.  & FWHM_res < 200. & FWHM_sky>0 & FWHM_res > 0);
  FWHM_res = FWHM_res(w);
  FWHM_sky = FWHM_sky(w);
  obsmode  = obsmode(w);
  m = max(max(FWHM_res),max(FWHM_sky));
  n = min(min(FWHM_res),min(FWHM_sky));
  winkill,1;window,1,dpi=90,style="aanda.gs";

  for(j=1;j<=numberof(mode);j++){
    w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, FWHM_res(w),FWHM_sky(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),.85*m,(1-j/10.)*m,color=c(j),tosys=1;
    }
  }
  
  plg, [n,m],[n,m],marks=0,type=2;
  xytitles,"Sky FWHM in H (mas)","Reconstructed FWHM (mas)";
  limits,n*.95,m*1.05;
  range,n*.95,m*1.05;
     
  if(is_string(aomode)){
    pltitle, meth + " in " +  aomode;
    pdf,"results/fwhmReconstructionIn" + aomode + "_" + method;
  }else{
    pltitle, meth;
    pdf,"results/statistics/fwhmReconstructionInAllAoModes" + "_" + method;
  }

  /*
  //////
  // Chi2
  ///////////////////////////
  winkill,2;window,2,dpi=90,style="aanda.gs";
  for(j=1;j<=numberof(mode);j++){
    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
      w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, chi2(w),0.103/r0(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),.85*m,(1-j/10.)*m,color=c(j),tosys=1;
    }
  }
  xytitles,"Seeing (arcsec)","!c^2";
  limits,0,1.5;range,0,nc*1.05;
  if(is_string(aomode)){
    pltitle, meth + " in " +  aomode;
    pdf,"results/chi2In" + aomode + "_" + method;
  }else{
    pltitle, meth;
    pdf,"results/statistics/chi2InAllAoModes" + "_" + method;
  }

  //////
  // Moffat
  ///////////////////////////
  winkill,3;window,3,dpi=90,style="aanda.gs";
  for(j=1;j<=numberof(mode);j++){
    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
      w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, I0res(w),I0sky(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),.85*m,(1-j/10.)*m,color=c(j),tosys=1;
    }
  }
  xytitles,"Sky I_0_","Reconstructed I_0_";
  range,0,.5;limits,0,.5;
  m = max(max(I0sky),max(I0res));
  plg, [0,m],[0,m],marks=0,type=2;
   
  winkill,4;window,4,dpi=90,style="aanda.gs";
  for(j=1;j<=numberof(mode);j++){
    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
      w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, ares(w),asky(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),.85*m,(1-j/10.)*m,color=c(j),tosys=1;
    }
  }
  xytitles,"Sky !a","Reconstructed !a";
  range,0,10.;limits,0,10.;
  m = max(max(asky),max(ares));
  plg, [0,m],[0,m],marks=0,type=2;

  
  winkill,5;window,5,dpi=90,style="aanda.gs";
  for(j=1;j<=numberof(mode);j++){
    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
      w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, bres(w),bsky(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),.85*m,(1-j/10.)*m,color=c(j),tosys=1;
    }
  }
  xytitles,"Sky !b","Reconstructed !b";
  range,0,4;limits,0.8,1.8;
  m = max(max(bsky),max(bres));
  plg, [0,m],[0,m],marks=0,type=2;
  */
}

func statisticsOnEE(method,aomode,n,all=)
/* DOCUMENT statisticsOnEE,"analytic",[],1.3,all=1
 */
{
  include,"pracConfig.i",1;
  
  savingDir  = "results/pracResults/resultsWith" + method + "_Vii_Method/"
  list  = listFile(savingDir);
  nfile = numberof(list);

  obsmode = EE_sky = EE_res = [];
  nPixelsCropped = 128;
  nPixels = 512;
  lambda = 1.677e-6;
  pixSize = 0.0297949;
  foV = nPixels * pixSize;
  D = 4.2;
  
  dl =  indgen(nPixelsCropped/2) * foV/nPixels;
  nn = where(dl>= n*radian2arcsec*lambda/D)(1);
 
  
  for(i=1;i<=nfile;i++){
    //loading results
    p = readfits(savingDir + list(i));
    //grabbing aomode
    mode = (strchar(*p(1)))(2);

    if(mode == aomode || all == 1){
      obsmode = grow(obsmode,mode);
      //grabbing the Strehl ratios
      EE_sky = grow(EE_sky,(*p(7))(nn,1));
      EE_res = grow(EE_res,(*p(7))(nn,2)) ;
      write,format="Files grabbed :%.3g%s\r",100.*i/nfile,"%";
    }
  }

  //mode = sortLabel(obsmode);
  mode = ["SCAO","MOAO4L3N","MOAO3N","MOAO4L3T","GLAO"];
  c  = ["black","red","blue","magenta","green"];
  
  mm = max(max(EE_res),max(EE_sky));
  mn = min(min(EE_res),min(EE_sky));
  winkill,0;window,0,dpi=90,style="aanda.gs";

  for(j=1;j<=numberof(mode);j++){

    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
      w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, EE_res(w),EE_sky(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),mn+mm/10.,(1-j/10.)*mm,color=c(j),tosys=1;
    }
  }

  plg, [mn,mm],[mn,mm],marks=0,type=2;
  xytitles,"Sky Ensquared energy (%) at "+var2str(n)+ "!l/D in H","Rec. Ensquared Energy (%) at "+var2str(n)+"!l/D";

  meth = "Estimation method";
  if(method == "ls")
    meth = "TS-based method";
  else if(method == "rtc")
    meth = "RTC-based method";

  nn = var2str(n);
  if(!is_integer(n)){
    nn = var2str(n);
    w  = strfind(".",nn);
    nn = streplace(nn,w,"-");
  }
  
  if(is_string(aomode)){
    pltitle, meth + " in " +  aomode;
    pdf,"results/eeReconstructionIn" + aomode + "_"+nn + "lonD_" + method;
  }else{
    pltitle, meth;
    pdf,"results/statistics/eeReconstructionInAllAoModes_"+nn + "lonD_" + method;
  }
}

/*
 ____ ___ ____  ____  _        _ __   _____ _   _  ____ 
|  _ \_ _/ ___||  _ \| |      / \\ \ / /_ _| \ | |/ ___|
| | | | |\___ \| |_) | |     / _ \\ V / | ||  \| | |  _ 
| |_| | | ___) |  __/| |___ / ___ \| |  | || |\  | |_| |
|____/___|____/|_|   |_____/_/   \_\_| |___|_| \_|\____|
                                                        

*/


func plotEE(timedata,procDir)
/* DOCUMENT plotEE,"02h31m11s","2013_09_17_onsky/"

 */
{
  p = pracMain(timedata,Dir=procDir,psfrMethod = "analytic",averageMode = "Vii",verb=0,disp=1);
  eesky = (*p(7))(,1);
  eea   = (*p(7))(,2);
  otfsky  = (*p(8))(,,1);
  otfa  = (*p(8))(,,2);
  window,1;pdf,"results/psf_analytic_" + timedata + ".pdf"; 
  p = pracMain(timedata,Dir=procDir,psfrMethod = "mmse",averageMode="",verb=0,disp=1);
  eei   = (*p(7))(,2);
  otfi  = (*p(8))(,,2);
  window,1;pdf,"results/psf_mmse_" + timedata + ".pdf"; 
  p = pracMain(timedata,Dir=procDir,psfrMethod = "ls",averageMode = "Vii",verb=0,disp=1);
  eev   = (*p(7))(,2);
  otfv  = (*p(8))(,,2);
  window,1;pdf,"results/psf_ls_" + timedata + ".pdf";
  
  boxsize = span(1,cam.nPixelsCropped-3.,cam.nPixelsCropped) * cam.pixSize;
  otel = roll(OTF_telescope(tel.diam,tel.obs,cam.nPixelsCropped,tel.pixSize*tel.nPixels/cam.nPixelsCropped));
  
  //Displaying EE
  winkill,9;window,9,style="aanda.gs",dpi=90;clr;

  plg, eesky, boxsize;
  plg, eea,boxsize,color=[64,64,64];
  plg, eev, boxsize,color=[100,100,100];
  plg, eei, boxsize,color=[128,128,128];
  plg, [100,100],[-0.1,max(boxsize)*1.05],type=2,marks=0;

  fcut = radian2arcsec*cam.lambda/tel.pitch;
  plg, [0,100],[fcut,fcut],type=2,marks=0;
  plg, [0,100],[1,1]*1.22*radian2arcsec*cam.lambda/tel.diam,type=2,marks=0;
  xytitles,"Angular separation from center (arcsec)","Ensquared Energy (%)";
  plt,"A: On-sky PSF",1.,40,tosys=1;
  plt,"B: Full analytic PSF",1.,30,tosys=1;
  plt,"C: TS-based LS PSF",1.,20,tosys=1;
  plt,"D: TS-based MMSE PSF",1.,10,tosys=1;
  plt,"DM cut frequency",fcut*1.05,50,tosys=1;
  plt,"1.22 !l/D",1.22*radian2arcsec*cam.lambda/tel.diam*1.05,10,tosys=1;
  range,0,105;
  limits,-0.1,max(boxsize)*1.05;
pdf,"results/ee_" + timedata + ".pdf";

  //Displaying OTFs
  winkill,10;window,10,style="aanda.gs",dpi=90;clr;
  dl =  indgen(cam.nPixelsCropped/2) * tel.foV/tel.nPixels;
  plg,otfsky(cam.nPixelsCropped/2+1,cam.nPixelsCropped/2+1:)/max(otfsky),dl;
  plg,otfa(cam.nPixelsCropped/2+1,cam.nPixelsCropped/2+1:)/max(otfa),dl;
  plg,otfi(cam.nPixelsCropped/2+1,cam.nPixelsCropped/2+1:)/max(otfi),dl;
  plg,otfv(cam.nPixelsCropped/2+1,cam.nPixelsCropped/2+1:)/max(otfv),dl;
  plg,otel(cam.nPixelsCropped/2+1,cam.nPixelsCropped/2+1:)/max(otel),dl,type=2,marks=0;

  xytitles,"D/!l","Normalized OTF";
  plt,"Dashed line: Perfect telescope",1,1,tosys=1;
  plt,"A: Sky OTF",1,0.8,tosys=1;
  plt,"B: Full analytic OTF",1,0.6,tosys=1;
  plt,"C: TS-based LS OTF",1,0.4,tosys=1;
  plt,"D:  TS-based MMSE OTF",1,0.2,tosys=1;
  range,-.1,1.1;
  pdf,"results/otfs_" + timedata + ".pdf";
}

func plotProfiles(void)
/* DOCUMENT plotProfiles;
  
 */
{

  cnh = readfits("results/Profiles/seeingProfile.fits");
  l0h = readfits("results/Profiles/l0hProfile.fits");
  vh = readfits("results/Profiles/vhProfile.fits");
  alt = readfits("results/Profiles/altProfile.fits");
  cnh0 = readfits("results/Profiles/seeing.fits");
  cnhg = readfits("results/Profiles/groundSeeing.fits");
  cnhalt = readfits("results/Profiles/altSeeing.fits");
  l0eff = readfits("results/Profiles/effL0.fits");
  l0g = readfits("results/Profiles/groundL0.fits");
  l0alt = readfits("results/Profiles/altL0.fits");
  v0 = readfits("results/Profiles/effV0.fits");
  v0g = readfits("results/Profiles/groundV0.fits");
  v0alt = readfits("results/Profiles/altV0.fits");
  
  if(0){
    //initialiazing vectors...
    cnh = alt = l0h = vh = dirh = [];//fitted values
    dcnh = dalt = dl0h = dvh = ddirh = [];//uncertainties

    listdir = listFile("profiles/");
    ndir = numberof(listdir);
    //managing directories
    for(j=1;j<=ndir;j++){
      dir      = listdir(j);
      pathdata = "profiles/" + dir + "/";
      listtl = listFile(pathdata);
      ntl    = numberof(listtl);
      dataDirRoot = dataDir + dir + "/";
      write,format="Processing directory %s\n",dir;
      //loop on file
      for(i=1;i<=ntl;i++){
        //determining the airmass
        timetl = strpart(listtl(i),21:29);
        restorefits,"slopestl",timetl,pathtl,fake=1;
        object = strcase(1,readFitsKey(pathtl,"OBJECT"));
        tmp = givesCoordinatesFromObject(object);
        array_ra = str2flt(decoupe(tmp(1),' '));
        array_dec = str2flt(decoupe(tmp(2),' '));
        airm = airmassFromDate(extractDate(pathtl),array_ra,array_dec);
        //growing results
        ptrprof = readfits(pathdata + listtl(i));
        cnh     = grow(cnh,*ptrprof(1)/airm);
        alt     = grow(alt,*ptrprof(2)/airm);
        l0h     = grow(l0h,*ptrprof(3));
        track   = grow(track,*ptrprof(4));
        vh      = grow(vh,*ptrprof(5));
        dirh    = grow(dirh,*ptrprof(6));
        dcnh    = grow(dcnh,*ptrprof(7)/airm);
        dalt    = grow(dalt,*ptrprof(8)/airm);
        dl0h    = grow(dl0h,*ptrprof(9));
        dtrack  = grow(dtrack,*ptrprof(10));
        dvh     = grow(dvh,*ptrprof(11));
        ddirh   = grow(ddirh,*ptrprof(12));
      
        write,format="Loading file:%.3g%s\r",100.*i/ntl,"%";
      
      }
    }
  
    
    alt    /= 1000;
    //computing integrated value
    cnh0 = l0eff = v0 = cnhalt = l0alt = v0alt = cnhg = l0g = v0g = [];
    for(i=1;i<=numberof(cnh)/5;i++){
      cnhi   = cnh(1+(i-1)*5:5*i);
      cnha   = cnhi(2:);
      l0hi   = l0h(1+(i-1)*5:5*i);
      wi     = where(l0hi<=100.)
        vhi    = vh(1+(i-1)*5:5*i);
      //global values
      cnh0   = grow(cnh0,0.103/sum(cnhi)^(-3/5.));
      if(is_array(wi))
        l0eff  = grow(l0eff,(sum(l0hi(wi)^(5/3.)*cnhi(wi))/sum(cnhi(wi)))^(3/5.));
      v0     = grow(v0,(sum(cnhi*vhi^(5/3.))/sum(cnhi))^(3/5.));
      //ground values
      cnhalt = grow(cnhalt,0.103/sum(cnha)^(-3/5.));
      if(numberof(wi) == 1 && wi(1) != 1){
        l0alt  = grow(l0alt,(sum(cnhi(wi)*l0hi(wi)^(5/3.))/sum(cnhi(wi)))^(3/5.));
      }else if(is_array(wi) && numberof(wi)!=1){
        l0alt  = grow(l0alt,(sum(cnhi(wi(2:))*l0hi(wi(2:))^(5/3.))/sum(cnhi(wi(2:))))^(3/5.));
      }
      v0alt  = grow(v0alt,(sum(cnha*vhi(2:)^(5/3.))/sum(cnha))^(3/5.));
      //values integrated onto altitude
      cnhg   = grow(cnhg,0.103/cnhi(1)^(-3/5.));
      l0g    = grow(l0g,l0hi(1));
      v0g    = grow(v0g,vhi(1));
    }

    cnh     = 0.103/cnh^(-3/5.);
    dcnh    = 0.103*3/5.*dcnh/cnh^(2/5.); 
    //wa      = where(alt>=20);
    //alt(wa) = 20.;
    w0 = where(l0h<=40);
    write,format="%.3g%s of jerk outer scale\n",100.*(1.- double(numberof(w0))/numberof(l0h)),"%s";

    writefits,"seeingProfile.fits",cnh;
    writefits,"l0hProfile.fits",l0h;
    writefits,"vhProfiles.fits",vh;
    writefits,"altProfile.fits",alt;
    writefits,"seeing.fits",cnh0;
    writefits,"groundSeeing.fits",cnhg;
    writefits,"altSeeing.fits",cnhalt;
    writefits,"effL0.fits",l0eff;
    writefits,"groundL0.fits",l0g;
    writefits,"altL0.fits",l0alt;
    writefits,"effV0.fits",v0;
    writefits,"groundV0.fits",v0g;
    writefits,"altV0.fits",v0alt;
  }
  write,format="%d processed files\n",numberof(cnh)/5;
  //..... DISPLAY .....//
  
  //Plotting altitude versus cnh
  winkill,0;window,0,dpi=90,style="aanda.gs";
  plmk,alt,cnh,marker=4,msize=0.1;
  range,-1,25;
  limits,0,2.;
  xytitles,"Seeing at 500 nm (arcsec)","Altitude (km)";
  pdf, "cnh_profile.pdf";

  //Plotting seeing histograms
  winkill,4;window,4,dpi=90,style="aanda.gs";
  res = .04;cnhlim = 2.;
  ns = int(minmax(cnh0(where(cnh0<cnhlim)))(dif)(1)/res-1.);
  histo,cnh0(where(cnh0<cnhlim)),ns;
  histo,cnhg(where(cnhg<cnhlim)),ns;
  histo,cnhalt(where(cnhalt<cnhlim)),int(ns/1.05);
  xytitles,"Seeing at 500 nm (arcsec)","Counts";
  lim = limits();
  plt,"A: Global seeing",lim(2)/2,lim(4)/1.5,tosys=1;
  plt,"B: Ground value ",lim(2)/2,lim(4)/1.5-lim(4)/8,tosys=1;
  plt,"C: Altitude value",lim(2)/2,lim(4)/1.5-lim(4)/4,tosys=1;
  limits,-.1,2.2;
  pdf, "cnh_histo.pdf";
  w0 = where(cnh0<=cnhlim);
  write,format="Average/1-sigma global seeing         = %.3g''/ %.3g''\n",avg(cnh0(w0)),(cnh0(w0))(rms);
  w0 = where(cnhg<=cnhlim);
  write,format="Average/1-sigma ground seeing         = %.3g''/ %.3g''\n",avg(cnhg(w0)),(cnhg(w0))(rms);
  w0 = where(cnhalt<=cnhlim);
  write,format="Average/1-sigma altitude seeing       = %.3g''/ %.3g''\n",avg(cnhalt(w0)),(cnhalt(w0))(rms);
  
  //Plotting altitude versus l0h
  winkill,1;window,1,dpi=90,style="aanda.gs";
  w0 = where(l0h<=100.);
  plmk,alt(w0),l0h(w0),marker=4,msize=0.1;
  range,-1,25;
  limits,-5,100;
  xytitles,"L_0_(h) (m)","Altitude (km)";
  pdf, "L0h_profile.pdf";

  //Plotting outer scale histograms
  winkill,5;window,5,dpi=90,style="aanda.gs";
  res = 1.;l0lim=100.;
  ns = int(minmax(l0eff(where(l0eff<=l0lim)))(dif)/res-1.);
  histo,l0eff(where(l0eff<=l0lim)),ns;
  histo,l0g(where(l0g<=l0lim)),ns;
  histo,l0alt(where(l0alt<=l0lim)),ns;
  xytitles,"Outer scale (m)","Counts";
  lim = limits();
  plt,"A: Effective outer scale",lim(2)/2,lim(4)/2,tosys=1;
  plt,"B: Ground value ",lim(2)/2,lim(4)/2-lim(4)/8,tosys=1;
  plt,"C: Altitude value",lim(2)/2,lim(4)/2-lim(4)/4,tosys=1;
  pdf, "l0h_histo.pdf";
  w0 = where(l0eff<=l0lim);
  write,format="Average/1-sigma effective outer scale = %.3g m /%.3g m\n",avg(l0eff(w0)),(l0eff(w0))(rms);
  w0 = where(l0g<=l0lim);
  write,format="Average/1-sigma ground outer scale    = %.3g m /%.3g m\n",avg(l0g(w0)),(l0g(w0))(rms);
  w0 = where(l0alt<=l0lim);
  write,format="Average/1-sigma altitude outer scale  = %.3g m /%.3g m\n",avg(l0alt(w0)),(l0alt(w0))(rms);

  //Plotting altitude versus vh
  vlim=14;
  winkill,2;window,2,dpi=90,style="aanda.gs";
  plmk,alt(where(vh<vlim)),vh(where(vh<vlim)),marker=4,msize=0.1;
  range,-1,25;
  limits,-1,vlim;
  xytitles,"v(h) (m/s)","Altitude (km)";
  pdf, "vh_profile.pdf";

  //Plotting wind speed histograms
  winkill,6;window,6,dpi=90,style="aanda.gs";
  res = .4;
  ns = int(minmax(v0(where(v0<vlim)))(dif)/res-1.);
  histo,v0(where(v0<vlim)),ns;
  histo,v0g(where(v0g<vlim)),ns;
  histo,v0alt(where(v0alt<vlim)),ns;
  xytitles,"Wind speed (m/s)","Counts";
  lim = limits();
  plt,"A: Temporal coherence wind velocity",lim(2)/2,lim(4)/2,tosys=1;
  plt,"B: Ground value ",lim(2)/2,lim(4)/2-lim(4)/8,tosys=1;
  plt,"C: Altitude value",lim(2)/2,lim(4)/2-lim(4)/4,tosys=1;
  pdf, "vh_histo.pdf";
  
  w0 = where(v0g<vlim);
  write,format="Average/1-sigma global wind speed     = %.3g m/s /%.3g m/s\n",avg(v0(w0)),(v0(w0))(rms);
  w0 = where(v0g<vlim);
  write,format="Average/1-sigma ground wind speed     = %.3g m/s /%.3g m/s\n",avg(v0g(w0)),(v0g(w0))(rms);
  w0 = where(v0alt<vlim);
  write,format="Average/1-sigma altitude wind speed   = %.3g m/s /%.3g m/s\n",avg(v0alt(w0)),(v0alt(w0))(rms);
}

/*
 _____ ____  ____   ___  ____  
| ____|  _ \|  _ \ / _ \|  _ \ 
|  _| | |_) | |_) | | | | |_) |
| |___|  _ <|  _ <| |_| |  _ < 
|_____|_| \_\_| \_\\___/|_| \_\
                               
 ____  ____  _____    _    _  ______   _____        ___   _ 
| __ )|  _ \| ____|  / \  | |/ /  _ \ / _ \ \      / / \ | |
|  _ \| |_) |  _|   / _ \ | ' /| | | | | | \ \ /\ / /|  \| |
| |_) |  _ <| |___ / ___ \| . \| |_| | |_| |\ V  V / | |\  |
|____/|_| \_\_____/_/   \_\_|\_\____/ \___/  \_/\_/  |_| \_|
                                                            
*/
func plotScriptErrorBreakdown(Dir,scriptsuffix,path_out,redo=,wb=)
/* DOCUMENT plotScriptErrorBreakdown,"2013_09_17_onsky/","script293","results/Scripts/";
 */
{
  local res;
  if(!redo)
    res = readfits(path_out + "scriptProcess_" + scriptsuffix +".fits",err=1);
  
  include, "pracConfig.i",1;
  dataDirRoot = dataDir + Dir;
  //takes all slopestl files
  pathstl = listVersion(dataDirRoot,"fits","slopestl");
  pathstl = pathstl(where(strfind(scriptsuffix,pathstl)(2,) != -1));
  ntl = numberof(pathstl);
    
  if(is_void(res)){
    s_ir = s_tomo = s_noise = s_bw = s_fit = s_static = s_ncpa = s_ol = [];
    aomode = srsky = srmar = srpar = srborn  = mode = [];
    r0 = L0 = v =  h = ved1 = ved2 = ved3 = ved4 = ved5 = [];
    cnh2 = alt2 = l0h2 = vh2 = [];
    PSF = array(0.,128,128,ntl);
    res = array(pointer,24);
    
    for(j=1;j<=ntl;j++){
      timetl  = extractDate(pathstl(j));
      h       = grow(h,str2time(strpart(timetl,12:)));
      tmp     = pracMain(timetl,Dir = Dir,budgetonly=1);
      mode    = grow(mode,strchar(*tmp(1))(2));
      
      PSF(,,j)  = (*tmp(6))(,,1);
      
      p       = *tmp(9);

      s_ir  = grow(s_ir,p(1));
      s_fit  = grow(s_fit,p(2));
      s_bw  = grow(s_bw,p(3));
      s_tomo  = grow(s_tomo,p(4));
      s_noise  = grow(s_noise,p(5));
      s_alias  = grow(s_alias,p(6));
      s_static = grow(s_static,p(7));
      s_ncpa = grow(s_ncpa,115);//grow(s_ncpa,p(8));
      s_ol  = grow(s_ol,p(9));
      
      p       = 100*(*tmp(11));
      srsky = grow(srsky,p(1));
      srmar = grow(srmar,p(3));
      srpar = grow(srpar,p(4));
      srborn = grow(srborn,p(5));
    
      r0 = grow(r0,(*tmp(2))(1));
      L0 = grow(L0,(*tmp(2))(2));
      v = grow(v,(*tmp(2))(3));

      ved1 = grow(ved1,(*tmp(10))(1));
      ved2 = grow(ved2,(*tmp(10))(2));
      ved3 = grow(ved3,(*tmp(10))(3));
      ved4 = grow(ved4,(*tmp(10))(4));
      ved5 = grow(ved5,(*tmp(10))(5));

      cnh2 = grow(cnh2,(*tmp(3))(,1));
      alt2 = grow(alt2,(*tmp(3))(,2));
      l0h2 = grow(l0h2,(*tmp(3))(,3));
      vh2  = grow(vh2,(*tmp(3))(,4));
      
    }
    res(1) = &strchar(mode);res(2) = &PSF;
    res(3) = &s_ir; res(4) = &s_tomo; res(5) = &s_alias; res(6) = &s_noise;
    res(7) = &s_bw;res(8) = &s_fit; res(9) = &s_ol;res(10) = &s_static;
    res(11) = &s_ncpa; res(12) = &srsky; res(13) = &srmar; res(14) = &srpar;res(15) = &srborn;
    res(16) = &r0;res(17) = &L0; res(18) = &v;res(19) = &h;
    res(20) = &[ved1,ved2,ved3,ved4,ved5];
    res(21) = &cnh2;
    res(22) = &alt2;
    res(23) = &l0h2;
    res(24) = &vh2;

    writefits,path_out + "scriptProcess_" + scriptsuffix +".fits",res;
  }

  ////////////
  //Displaying error breakdown
  /////////////////////////////////////////////////

  if(wb!=1){
    c1 = "blue";
    c2 = "black";
    c3 = "red";
    col1=[char(249)];
    col2=[char(241)];
    col3=[char(251)];
  }else{
    c1 = [96,96,96];
    c2 = "black";
    c3 = [128,128,128];
    col1=[char(242)];
    col2=[char(241)];
    col3=[char(243)];
  }
  savedir = "/home/omartin/Documents/Articles/Canary_Phase_B/Figures/" + scriptsuffix +"/";
  if(!direxist(savedir))
    system,"mkdir " + savedir;

  
  modes = strchar(*res(1));
  sir_1  = (*res(3))(1::3);sir_2  = (*res(3))(2::3);sir_3  = (*res(3))(3::3);
  stomo_1  = (*res(4))(1::3);stomo_2  = (*res(4))(2::3);stomo_3  = (*res(4))(3::3);
  salias_1  = (*res(5))(1::3);salias_2  = (*res(5))(2::3);salias_3  = (*res(5))(3::3);
  snoise_1  = (*res(6))(1::3);snoise_2  = (*res(6))(2::3);snoise_3  = (*res(6))(3::3);
  sbw_1  = (*res(7))(1::3);sbw_2  = (*res(7))(2::3);sbw_3  = (*res(7))(3::3);
  sfit_1  = (*res(8))(1::3);sfit_2  = (*res(8))(2::3);sfit_3  = (*res(8))(3::3);
  sol_1  = (*res(9))(1::3);sol_2  = (*res(9))(2::3);sol_3  = (*res(9))(3::3);
  sstatic_1  = (*res(10))(1::3);sstatic_2  = (*res(10))(2::3);sstatic_3  = (*res(10))(3::3);
  sncpa_1  = (*res(11))(1::3);sncpa_2  = (*res(11))(2::3);sncpa_3  = (*res(11))(3::3);
  
  y1 = [sqrt(avg(sir_1^2)),sqrt(avg(stomo_1^2)),sqrt(avg(salias_1^2)),sqrt(avg(snoise_1^2)),sqrt(avg(sbw_1^2)),sqrt(avg(sfit_1^2)),sqrt(avg(sol_1^2)),sqrt(avg(sstatic_1^2)),sqrt(avg(sncpa_1^2))];
  y2 = [sqrt(avg(sir_2^2)),sqrt(avg(stomo_2^2)),sqrt(avg(salias_2^2)),sqrt(avg(snoise_2^2)),sqrt(avg(sbw_2^2)),sqrt(avg(sfit_2^2)),sqrt(avg(sol_2^2)),sqrt(avg(sstatic_2^2)),sqrt(avg(sncpa_2^2))];
  y3 = [sqrt(avg(sir_3^2)),sqrt(avg(stomo_3^2)),sqrt(avg(salias_3^2)),sqrt(avg(snoise_3^2)),sqrt(avg(sbw_3^2)),sqrt(avg(sfit_3^2)),sqrt(avg(sol_3^2)),sqrt(avg(sstatic_3^2)),sqrt(avg(sncpa_3^2))];
  
  labs = ["!s_!e-IT","!s_Tomography","!s_Aliasing","!s_Noise","!s_ServoLag","!s_Fitting","!s_Go-to","!s_Static","!s_ncpa"];
  winkill,0;window,0,style="aanda.gs",dpi=90;clr;gridxy,1,1;
  step=2;
  plotsBarDiagram,y1,y2=y2,y3=y3,labs,col1=col1,col2=col2,col3=col3,title=0,step=step;
  plt,modes(1),7*step,.8*max(y1(1),y2(1),y3(1)),tosys=1,color=c1;
  plt,modes(2),7*step,.7*max(y1(1),y2(1),y3(1)),tosys=1,color=c2;
  plt,modes(3),7*step,.6*max(y1(1),y2(1),y3(1)),tosys=1,color=c3;
  xytitles," ","Wave Front Error [nm rms]";
  pdf,savedir + scriptsuffix + "_budget.pdf";

  /////////
  //Displaying averaged PSF
  /////////////////////////////////////////////////

  
  PSF1 = ((*res(2))(,,1::3))(,,avg);
  PSF2 = ((*res(2))(,,2::3))(,,avg);
  PSF3 = ((*res(2))(,,3::3))(,,avg);
  a=3. ;
  cutmax = (max(max(PSF1),max(PSF2),max(PSF3)))^(1./a);
  box = dimsof(PSF1)(0); uz = 0.0297949; 

  window,1;clr;//palette,"gray.gp";
  pli,(abs(PSF1))^(1./a),-box*uz/2,-box*uz/2,box*uz/2,box*uz/2,cmin = 0,cmax = cutmax;
  xytitles,"Arcsec","Arcsec";
  pltitle,"On-sky PSF in " + modes(1) + "\n Strehl = " + var2str(arrondi(100*max(PSF1),1)) + "%";
  pdf,savedir+ scriptsuffix +"_" + modes(1)+ "_psf.pdf";

  window,2;clr;//palette,"gray.gp";
  pli,(abs(PSF2))^(1./a),-box*uz/2,-box*uz/2,box*uz/2,box*uz/2,cmin = 0,cmax = cutmax;
  xytitles,"Arcsec","Arcsec";
  pltitle,"On-sky PSF in " + modes(2) + "\n Strehl = " + var2str(arrondi(100*max(PSF2),1)) + "%";
  pdf,savedir + scriptsuffix +"_" + modes(2)+ "_psf.pdf";

  window,3;clr;//palette,"gray.gp";
  pli,(abs(PSF3))^(1./a),-box*uz/2,-box*uz/2,box*uz/2,box*uz/2,cmin = 0,cmax = cutmax;
  xytitles,"Arcsec","Arcsec";
  pltitle,"On-sky PSF in " + modes(3) + "\n Strehl = " + var2str(arrondi(100*max(PSF3),1)) + " %";
  pdf,savedir+ scriptsuffix +"_" + modes(3)+ "_psf.pdf";

  /////////
  //Displaying circular average of PSF/OTF
  /////////////////////////////////////////////////
  /*
  winkill,7;window,7,dpi=90,style="aanda.gs";clr;gridxy,1,1;//palette,"gray.gp";
  //plg,log(circularAveragePsf(PSFtel)),1e3*(indgen(box/2)-1)*uz,type=2;
  plg,log(abs(circularAveragePsf(PSF1))),1e3*(indgen(box/2)-1)*uz,color=c1;
  plg,log(abs(circularAveragePsf(PSF2))),1e3*(indgen(box/2)-1)*uz,color=c2;
  plg,log(abs(circularAveragePsf(PSF3))),1e3*(indgen(box/2)-1)*uz,color=c3;
  //plt,"A: Perfect telescope PSF",600,-1,tosys=1;
  m = log(max(max(PSF1),max(PSF2),max(PSF3)))
    plt,"A: " + modes(1) +" PSF",100,-8,tosys=1,color=c1;
  plt,"B: " +modes(2) + " PSF",100,-10,tosys=1,color=c2;
  plt,"C: " +modes(3) + " PSF",100,-12,tosys=1,color=c3;
  plt,scriptsuffix,100,-14,tosys=1;
  xytitles,"Angular separation [mas]","Normalized PSF (log scale)";
  limits,-10,1000;
  pdf,savedir + scriptsuffix +"_circPSF.pdf";

  winkill,8;window,8,dpi=90,style="aanda.gs";clr;gridxy,1,1;//palette,"gray.gp";
  OTF1 = roll(fft(roll(PSF1)).re);
  OTF2 = roll(fft(roll(PSF2)).re);
  OTF3 = roll(fft(roll(PSF3)).re);
  dl =  ((indgen(box/2)-1) * uz*2);
  plg,log(abs(circularAveragePsf(OTF1/max(OTF1)))),dl,color=c1;
  plg,log(abs(circularAveragePsf(OTF2/max(OTF2)))),dl,color=c2;
  plg,log(abs(circularAveragePsf(OTF3/max(OTF3)))),dl,color=c3;
  //plt,"A: Perfect telescope PSF",600,-1,tosys=1;
  plt,"A: " +modes(1) +" OTF",.1,-4,tosys=1,color=c1;
  plt,"B: " +modes(2) + " OTF",.1,-5,tosys=1,color=c2;
  plt,"C: " +modes(3) + " OTF",.1,-6,tosys=1,color=c3;
  plt,scriptsuffix,.1,-7,tosys=1;
  xytitles,"Normalized Frequency [!r/!l units]","Normalized OTF (log scale)";
  limits,-.1,1.3;range,-10,0;
  pdf,savedir + scriptsuffix +"_circOTF.pdf";
  */
  
  ////////////
  //Displaying Strehl ratios versus time
  /////////////////////////////////////////////////

  w1 = where(modes == modes(1))
  w2 = where(modes == modes(2));
  w3 = where(modes == modes(3));
  SRsky1 = (*res(12))(w1); SRsky2 = (*res(12))(w2); SRsky3 = (*res(12))(w3);
  SRmar1 = (*res(13))(w1); SRmar2 = (*res(13))(w2); SRmar3 = (*res(13))(w3);
  SRpar1 = (*res(14))(w1); SRpar2 = (*res(14))(w2); SRpar3 = (*res(14))(w3);
  SRborn1 = (*res(15))(w1); SRborn2 = (*res(15))(w2); SRborn3 = (*res(15))(w3);
  h1 = (*res(19))(w1); h2 = (*res(19))(w2); h3 = (*res(19))(w3);

  if(SRsky1(avg) >= 50){
    SRanal1 = SRmar1(sort(h1));
  }else if (50 > SRsky1(avg)  && SRsky1(avg) >= 15){
    SRanal1 = SRpar1(sort(h1));
  }else if (15 > SRsky1(avg)){
    SRanal1 = SRborn1(sort(h1));
  }

  if(SRsky2(avg) >= 50){
    SRanal2 = SRmar2(sort(h2));
  }else if (50 > SRsky2(avg) && SRsky2(avg) >= 15){
    SRanal2 = SRpar2(sort(h2));
  }else if (15 > SRsky2(avg)){
    SRanal2 = SRborn2(sort(h2));
  }

  if(SRsky3(avg) >= 50){
    SRanal3 = SRmar3(sort(h3));
  }else if (50 > SRsky3(avg) && SRsky3(avg) >= 15){
    SRanal3 = SRpar3(sort(h3));
  }else if (15 > SRsky3(avg)){
    SRanal3 = SRborn3(sort(h3));
  }

  SRanal1 = SRpar1(sort(h1));SRanal2 = SRpar2(sort(h2));SRanal3 = SRpar3(sort(h3));
  
  winkill,4;window,4,style="aanda.gs",dpi=90;clr;gridxy,1,1;
  plmk, SRsky1(sort(h1)),h1(sort(h1)),msize=.4,marker=3,width=10,color=c1;
  plg,SRsky1(sort(h1)),h1(sort(h1)),marks=0,color=c1;
  plmk, SRanal1,h1(sort(h1)),msize=.4,marker=3,color=c1;
  plg,SRanal1,h1(sort(h1)),marks=0,type=2,color=c1;

  plmk, SRsky2(sort(h2)),h2(sort(h2)),msize=.4,marker=4,width=10,color=c2;
  plg,SRsky2(sort(h2)),h2(sort(h2)),marks=0,color=c2;
  plmk, SRanal2,h2(sort(h2)),msize=.4,marker=4,color=c2;
  plg,SRanal2,h2(sort(h2)),marks=0,type=2,color=c2;

  plmk, SRsky3(sort(h3)),h3(sort(h3)),msize=.4,marker=6,width=10,color=c3;
  plg,  SRsky3(sort(h3)),h3(sort(h3)),marks=0,color=c3;
  plmk, SRanal3,h3(sort(h3)),msize=.4,marker=6,color=c3;
  plg,  SRanal3,h3(sort(h3)),marks=0,type=2,color=c3;
  
  range,0,50;
  xytitles,"Local time [hour]","H-band SR [%]";
  mm = min(grow(h1,h2,h3));
  mn = max(grow(h1,h2,h3));
  plt,"Triangles: " + modes(1),mm,48,tosys=1,color=c1;
  plt,"Circles: " + modes(2),mm,45,tosys=1,color=c2;
  plt,"Crosses: " + modes(3),mm,42,tosys=1,color=c3;
  plt,"Solid line: on-sky SR",mn-.05,48,tosys=1;
  plt,"Dash  line: analytic SR",mn-.05,45,tosys=1;
  plt,"Script " + strpart(scriptsuffix,7:),mn-.05,42,tosys=1;
  limits,mm-0.01,mn+0.01;
  pdf,savedir + scriptsuffix +"_srvshour.pdf";

  //////////
  //Displaying VED
  /////////////////////////////////////////////////

   cnh = [median((*res(21))(1::5)),median((*res(21))(2::5)),median((*res(21))(3::5)),median((*res(21))(4::5)),median((*res(21))(5::5))];
   alt = [median((*res(22))(1::5)),median((*res(22))(2::5)),median((*res(22))(3::5)),median((*res(22))(4::5)),median((*res(22))(5::5))]/1000.;
    
  w1 = where(modes!="SCAO");
  if(is_array(w1)){
    m = modes(w1);
    w2 = w1(where(m==m(1)));
    w3 = w1(where(m==m(2)));

    ved1 = sqrt(((*res(20))(w2,)^2)(avg,));
    ved2 = sqrt(((*res(20))(w3,)^2)(avg,));
   
  
    winkill,5;window,5,style="aanda.gs",dpi=90;clr;gridxy,1,1;
    thick=.1;
    for(j=1; j<=5; j++){
      d = 0.5*(0.5-thick);
      d = 2*thick;
      plfp, col2,[0,1,1,0]*ved1(j),([1,1,1,1]*alt(j) -d/2) + [-1,-1,1,1]*thick,[4];
      plfp, col3,[0,1,1,0]*ved2(j),([1,1,1,1]*alt(j) +d/2) + [-1,-1,1,1]*thick,[4];
    }

    xytitles,"Altitude [km]","VED [nm]";
    ml = max(grow(ved1,ved2));
    plt,m(1),1,ml*0.9,tosys=1,color=c2;
    plt,m(2),1,ml*0.8,tosys=1,color=c3;
    limits,min(alt)-1,max(max(alt)+1,15);
    range,0,ml*1.05;
    pdf,savedir + scriptsuffix + "_VED.pdf";
  
    ////////////
    // Displaying profiles
    /////////////////////////////////////////////////
    
    winkill,6;window,6,style="aanda.gs",dpi=90;clr;gridxy,1,1;
    wmoao = where(modes == "MOAO4L3N");
    suffcnh = readFitsKey(pathstl(wmoao(1)),"CN2H");
    ptrcnh  = restorefits("tomoparam",suffcnh);
    displayLayers,cnh,alt*1000,col=col3,thick=.2,percent=1;
    displayLayers,-(*ptrcnh(2)),*ptrcnh(3),col=[char(241)],percent=1,thick=.2;
    limits,-100,100;
    range,-1,22;
    plt,"Used profile\n to compute R",-90,15,tosys=1;
    plt,"Average profile\n during the script",20,15,tosys=1,color=[128,128,128];
    plt,"Script " + strpart(scriptsuffix,7:),-10,20,tosys=1;
    pdf,savedir + scriptsuffix + "_profiles.pdf";
  }
  //gives statistics on turbulence
  //r0
  r0 = avg((*res(16))^(-5/3.))^(-3/5.);
  cnh2 = cnh*r0^(-5/3.)/sum(cnh);
  r0g = sum(cnh2(where(alt <= 1)))^(-3/5.);
  r0a = sum(cnh2(where(alt >1)))^(-3/5.);
  //seeing
  seeing = 0.103/r0;
  sg    = 0.103/r0g;
  sa    = 0.103/r0a;
  //Wind speed
  vh = [median((*res(24))(1::5)),median((*res(24))(2::5)),median((*res(24))(3::5)),median((*res(24))(4::5)),median((*res(24))(5::5))];
  v  = avg(vh); 
  vg = avg(vh(where(alt<=0.1)));
  va = avg(vh(where(alt>0.1)));

  //Outerscale
  alt = (*res(22));
  l0h = (*res(23));
  w = where(l0h <=100.);
  l0h = l0h(w);
  alt = alt(w);
  L0 =  avg(l0h^(5/3.))^(3/5.);
  L0g = avg((l0h(where(alt<=0.1)))^(5/3.))^(3/5.);
  L0a = avg((l0h(where(alt>0.1)))^(5/3.))^(3/5.);

  write,format="r0 (500nm)       = %.3g cm\n",100*r0;
  write,format="Ground r0        = %.3g cm\n",100*r0g;
  write,format="Altitude r0      = %.3g cm\n",100*r0a;
  write,format="Seeing           = %.3g ''\n",seeing;
  write,format="Ground seeing    = %.3g ''\n",sg;
  write,format="Altitude seeing  = %.3g ''\n",sa;
  write,format="Wind speed       = %.3g m/s\n",v;
  write,format="Ground windspeed = %.3g m/s\n",vg;
  write,format="Alt. windspeed   = %.3g m/s\n",va;
  write,format="L0               = %.3g m\n",L0;
  write,format="Ground L0        = %.3g m\n",L0g;
  write,format="Altitude L0      = %.3g m\n",L0a;
  
}




  

/*
 ____  _  ____   __
/ ___|| |/ /\ \ / /
\___ \| ' /  \ V / 
 ___) | . \   | |  
|____/|_|\_\  |_|  
                   
 ____  _____ ____  _____ ___  ____  __  __    _    _   _  ____ _____ 
|  _ \| ____|  _ \|  ___/ _ \|  _ \|  \/  |  / \  | \ | |/ ___| ____|
| |_) |  _| | |_) | |_ | | | | |_) | |\/| | / _ \ |  \| | |   |  _|  
|  __/| |___|  _ <|  _|| |_| |  _ <| |  | |/ ___ \| |\  | |___| |___ 
|_|   |_____|_| \_\_|   \___/|_| \_\_|  |_/_/   \_\_| \_|\____|_____|
                                                                     
*/


func skyPerformance(aomode,all=)
/* DOCUMENT  skyPerformance,[],all=1;
 */
{

  include,"pracConfig.i",1;

  savingDir  = "results/pracResults/resultsWithls_Vii_Method/"
  list       = listFile(savingDir);
  nfile      = numberof(list);
  pixSize    =  0.0297949;
  nPixelsCropped = 128;
  nPixels = 512;
  lambda = 1.677e-6; 
  foV        = nPixels * pixSize;
  D          = 4.2;
  dl =  indgen(nPixelsCropped/2) * foV/nPixels;
  nn = where(dl>= 10*radian2arcsec*lambda/D)(1);
  
  obsmode = EE = SR = FWHM = seeing = [];
  for(i=1;i<=nfile;i++){
    //loading results
    p = readfits(savingDir + list(i));
    //grabbing aomode
    mode = strchar(*p(1))(2);

    if(mode == aomode || all == 1){
      obsmode = grow(obsmode,mode);
      //grabbing the Strehl ratios
      SR   = grow(SR,(*p(11))(1));
      ps   = (*p(6))(,,1);
      FWHM = grow(FWHM,1e3*getPsfFwhm(ps,pixSize,fit=0));
      EE   = grow(EE,(*p(7))(nn,1));
      seeing = grow(seeing, 0.103/(*p(2))(1));
      write,format="Files grabbed :%.3g%s\r",100.*i/nfile,"%";
    }
  }

  //mode = sortLabel(obsmode);
  mode = ["SCAO","MOAO4L3N","MOAO3N","MOAO4L3T","GLAO"];
  c  = ["black","red","blue","magenta","green"];
  

  winkill,0;window,0,dpi=90,style="aanda.gs";
  m =max(SR);
  m2 =max(seeing);
  for(j=1;j<=numberof(mode);j++){
    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
      w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, SR(w),seeing(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),1.5,(1-j/10.)*m,color=c(j),tosys=1;
    }
  }
  pdf,"results/statistics/skySR.pdf";
  
  xytitles,"Seeing at 500 nm(arcsec)", "Sky Strehl ratio in H";
  limits,0,2;
  range,-1,m*1.05;

  winkill,1;window,1,dpi=90,style="aanda.gs";
  m =max(FWHM);
  for(j=1;j<=numberof(mode);j++){
    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
      w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, FWHM(w),seeing(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),.1,(1-j/10.)*200,color=c(j),tosys=1;
    }
  }
  
  xytitles,"Seeing at 500 nm(arcsec)", "PSF FWHM in H";
  pdf,"results/statistics/skyFWHM.pdf";
  limits,0,2;
  range,70,200;
  
   winkill,2;window,2,dpi=90,style="aanda.gs";
  m =max(EE);
  for(j=1;j<=numberof(mode);j++){
    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
      w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, EE(w),seeing(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),.1,(1-j/10.)*m,color=c(j),tosys=1;
    }
  }
  
  xytitles,"Seeing at 500 nm(arcsec)","EE at 10!l/D (%)";
  pdf,"results/statistics/skyEE.pdf";
  limits,0,2;
  range,0,100;
}


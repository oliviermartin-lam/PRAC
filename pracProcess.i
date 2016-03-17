include, "pracMain.i";

//00h14
// r0              = 18.4 cm
// L0              = 9.45 m
// v               = 8.84 m/s
// Sky  Strehl     = 29.9 %
// TS   Strehl     = 29.6 %
// Est. Strehl     = 31 %
// RTC  Strehl     = 31.6 %
// Sky  FWHM       = 90.17 mas
// TS   FWHM       = 88.75 mas
// Est. FWHM       = 90.44 mas
// RTC  FWHM       = 90.32 mas
// TS   Chi2       = 3.295
// Est. Chi2       = 3.252
// RTC  Chi2       = 3.446
// Sky  Mof. prof. = [0.321, 1.73, 1.22]
// TS   Mof. prof. = [0.336, 1.82, 1.35]
// Est. Mof. prof. = [0.352, 2.09, 1.64]
// RTC  Mof. prof. = [0.358, 2.16, 1.73]

//00h15
// r0              = 20.1 cm
// L0              = 8.99 m
// v               = 6.12 m/s
// Sky  Strehl     = 24 %
// TS   Strehl     = 22.5 %
// Est. Strehl     = 22.4 %
// RTC  Strehl     = 21.1 %
// Sky  FWHM       = 106.6 mas
// TS   FWHM       = 100.3 mas
// Est. FWHM       = 98.54 mas
// RTC  FWHM       = 101.1 mas
// TS   Chi2       = 1.824
// Est. Chi2       = 1.74
// RTC  Chi2       = 1.808
// Sky  Mof. prof. = [0.247, 2.06, 1.24]
// TS   Mof. prof. = [0.246, 1.85, 1.15]
// Est. Mof. prof. = [0.247, 1.81, 1.14]
// RTC  Mof. prof. = [0.23, 1.82, 1.11]

//00h16
// r0              = 16.3 cm
// L0              = 8.14 m
// v               = 7.26 m/s
// Sky  Strehl     = 16.7 %
// TS   Strehl     = 15.7 %
// Est. Strehl     = 14.4 %
// RTC  Strehl     = 15.6 %
// Sky  FWHM       = 128.1 mas
// TS   FWHM       = 115.3 mas
// Est. FWHM       = 117.3 mas
// RTC  FWHM       = 112.4 mas
// TS   Chi2       = 1.344
// Est. Chi2       = 1.374
// RTC  Chi2       = 1.351
// Sky  Mof. prof. = [0.17, 2.43, 1.2]
// TS   Mof. prof. = [0.164, 1.99, 1.04]
// Est. Mof. prof. = [0.151, 1.99, 1.01]
// RTC  Mof. prof. = [0.165, 1.92, 1.03]
                   
  
func comparePSFRmethods(timedata,Dir=,disp=)
/* DOCUMENT comparePSFRmethods,timedata,Dir="2013_09_13_onsky/"

 */
{
  res_ts  = pracMain(timedata,Dir=Dir,psfrMethod="ts",disp=disp);
  res_est = pracMain(timedata,Dir=Dir,psfrMethod="estimation",disp=disp);
  res_rtc = pracMain(timedata,Dir=Dir,psfrMethod="rtc",disp=disp);

  EEsky  = (*res_ts(7))(,1);
  EEts   = (*res_ts(7))(,2);
  EEest  = (*res_est(7))(,2);
  EErtc  = (*res_rtc(7))(,2);

 
  // Relative Error on EE
  
  winkill,10;window,10,dpi=90,style="aanda.gs";
  logxy,1,0;gridxy,1,1;
  boxwidth = 1e3*span(1,cam.nPixelsCropped-3.,cam.nPixelsCropped) * cam.pixSize;
  plg, 100*(EEts-EEsky)/EEsky,boxwidth;
  plg, 100*(EEest-EEsky)/EEsky,boxwidth,color = [128,128,128];
  plg, 100*(EErtc-EEsky)/EEsky,boxwidth,color = [100,100,100];
  xytitles,"Integrated box width (mas)","Total EE error (%)";
  plt,"A: TS-based  PSF",.5e3,-5,tosys=1;
  plt,"B: Estimated PSF",.5e3,-6,tosys=1,color = [128,128,128];
  plt,"C: RTC-based PSF",.5e3,-7,tosys=1,color = [100,100,100];
  pltitle,rtc.aoMode + " at " + timedata;
  pdf,"results/EEerror_" + timedata;


  // Absolute Error on OTF
  
  os   = (*res_ts(8))(,,1);
  ots  = (*res_ts(8))(,,2);
  oest = (*res_est(8))(,,2);
  ortc = (*res_rtc(8))(,,2);
   
  os_avg   = circularAveragePsf(os);os_avg/=max(os_avg);
  ots_avg  = circularAveragePsf(ots);ots_avg/=max(ots_avg);
  oest_avg = circularAveragePsf(oest);oest_avg/=max(oest_avg);
  ortc_avg = circularAveragePsf(ortc);ortc_avg/=max(ortc_avg);

  ots_err  = ots_avg  - os_avg;
  oest_err = oest_avg - os_avg;
  ortc_err = ortc_avg - os_avg;
  
  winkill,11;window,11,style="aanda.gs",dpi=90;clr;
  dl =  (indgen(cam.nPixelsCropped/2) * tel.foV/tel.nPixels);
  n = numberof(dl);gridxy,1,1;
  plg, ots_err(1:n/2+5),dl(1:n/2+5);
  plg, oest_err(1:n/2+5),dl(1:n/2+5);
  plg, ortc_err(1:n/2+5),dl(1:n/2+5);
  xytitles,"!r/D","Absolute OTF error ";
  plt,"A: TS-based  OTF",.5,-.03,tosys=1;
  plt,"B: Estimated OTF",.5,-.04,tosys=1,color = [128,128,128];
  plt,"C: RTC-based OTF",.5,-.05,tosys=1,color = [100,100,100];
  pltitle,rtc.aoMode + " at " + timedata;
  
  pdf,"results/OTFerror_" + timedata;
  
  // Scalar results on performance

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

  write,format="Sky  Mof. prof. = [%.3g, %.3g, %.3g]\n",(*res_ts(12))(1,1),(*res_ts(12))(2,1),(*res_ts(12))(3,1);
  write,format="TS   Mof. prof. = [%.3g, %.3g, %.3g]\n", (*res_ts(12))(1,2),(*res_ts(12))(2,2),(*res_ts(12))(3,2);
  write,format="Est. Mof. prof. = [%.3g, %.3g, %.3g]\n",(*res_est(12))(1,2),(*res_est(12))(2,2),(*res_est(12))(3,2);
  write,format="RTC  Mof. prof. = [%.3g, %.3g, %.3g]\n", (*res_rtc(12))(1,2),(*res_rtc(12))(2,2),(*res_rtc(12))(3,2);

  // Moffat profile
  winkill,12;window,12,style="aanda.gs",dpi=90;clr;
  gridxy,1,1;
  y1  = (*res_ts(12))(,1);
  y2  = (*res_ts(12))(,2);
  y3  = (*res_est(12))(,2);
  y4  = (*res_rtc(12))(,2);
  tab = [y1,y2,y3,y4];
  multiBarDiagram,tab,["I_0","!a","!b"];
  xytitles," ","Moffat parameters";
  plt,"Sky PSF",.2,max(tab),tosys=1,color=[192,192,192];
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
  Dir = ["2013_09_13_onsky","2013_09_16_onsky","2013_09_15_onsky","2013_09_18_onsky", "2013_09_17_onsky"];

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
          
          if(method == "rtc" && mode == "SCAO"){
            continue;
          }else{
            p = pracMain(timetl,Dir=procDir,verb=0,disp=0,psfrMethod=method,writeRes=1);
          }

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

func statisticsOnSR(method,aomode,all=)
/* DOCUMENT statisticsOnSR,"analytic",[],all=1
 */
{

  include,"pracConfig.i",1;

  savingDir  = "results/pracResults/resultsWith" + method + "_Vii_Method/"
  list       = listFile(savingDir);
  nfile      = numberof(list);
  
  obsmode = SR_sky = SR_res = [];
  for(i=1;i<=nfile;i++){
    //loading results
    p = readfits(savingDir + list(i));
    //grabbing aomode
    mode = strchar(*p(1))(2);

    if(mode == aomode || all == 1){
      obsmode = grow(obsmode,mode);
      //grabbing the Strehl ratios
      SR_sky = grow(SR_sky,(*p(11))(1));
      SR_res = grow(SR_res,(*p(11))(2)) ;
      write,format="Files grabbed :%.3g%s\r",100.*i/nfile,"%";
    }
  }

  //mode = sortLabel(obsmode);
  mode = ["SCAO","MOAO4L3N","MOAO3N","MOAO4L3T","GLAO"];
  c  = ["black","red","blue","magenta","green"];
  m = max(max(SR_res),max(SR_sky));
  winkill,0;window,0,dpi=90,style="aanda.gs";

  for(j=1;j<=numberof(mode);j++){
    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
      w = where(strpart(obsmode,1:4) == mode(j));
    
    if(is_array(w)){
      plmk, SR_res(w),SR_sky(w),marker = j,msize=.2,width = 10,color=c(j);
      plt,mode(j),5,(1-j/10.)*m,color=c(j),tosys=1;
    }
  }
  
  plg, [0,m],[0,m],marks=0,type=2;
  xytitles,"Sky Strehl ratio in H","Reconstructed Strehl ratio";
  limits,-1,m*1.05;
  range,-1,m*1.05;
  
 meth = "Estimation method";
  if(method == "ls")
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
}

func statisticsOnFWHM(method,aomode,all=)
/* DOCUMENT statisticsOnFWHM,"analytic",[],all=1
 */
{

  include,"pracConfig.i",1;

  savingDir  = "results/pracResults/resultsWith" + method + "_Vii_Method/"
  list       = listFile(savingDir);
  nfile      = numberof(list);
  pixSize    = 0.0297949;
  //pixSize   *= 128./512.;
  
  obsmode = FWHM_sky = FWHM_res = [];
  os = or = array(0.,512,512);
  N  = 128;
  xi = 256-N/2+1;
  xf = 256 + N/2;
    
  for(i=1;i<=nfile;i++){
    //loading results
    p = readfits(savingDir + list(i));
    //grabbing aomode
    mode = strchar(*p(1))(2);

    if(mode == aomode || all == 1){
      obsmode = grow(obsmode,mode);
      //FWHM_sky = grow(FWHM_sky,(*p(11))(14));
      //FWHM_res = grow(FWHM_res,(*p(11))(13)) ;
      //grabbing PSFs
      ps = (*p(6))(,,1);
      pr = (*p(6))(,,2);

      /*
      //grabbing OTFs
      otfSky = (*p(8))(,,1);
      otfRes = (*p(8))(,,2);
      os(xi:xf,xi:xf) = otfSky;
      or(xi:xf,xi:xf) = otfRes;
      ps = roll(fft(roll(os)).re);
      pr = roll(fft(roll(or)).re);
      ps/=max(ps);pr/=max(pr);v
      */
      //computing FWHM
      FWHM_sky = grow(FWHM_sky,1e3*getPsfFwhm(ps,pixSize,fit=0));
      FWHM_res = grow(FWHM_res,1e3*getPsfFwhm(pr,pixSize,fit=0)) ;
      
      write,format="Files grabbed :%.3g%s\r",100.*i/nfile,"%";
    }
  }

  //mode = sortLabel(obsmode);
  mode = ["SCAO","MOAO4L3N","MOAO3N","MOAO4L3T","GLAO"];
  c  = ["black","red","blue","magenta","green"];
  w = where(FWHM_sky < 120.  & FWHM_res < 120. & FWHM_sky>0 & FWHM_res > 0);
  FWHM_res = FWHM_res(w);
  FWHM_sky = FWHM_sky(w);
  obsmode  = obsmode(w);
  m = max(max(FWHM_res),max(FWHM_sky));
  n = min(min(FWHM_res),min(FWHM_sky));
  winkill,0;window,0,dpi=90,style="aanda.gs";

  for(j=1;j<=numberof(mode);j++){
    if(mode(j) != "GLAO")
      w = where(obsmode == mode(j));
    else
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
  
 meth = "Estimation method";
  if(method == "ls")
    meth = "TS-based method";
  else if(method == "rtc")
    meth = "RTC-based method";
    
  if(is_string(aomode)){
    pltitle, meth + " in " +  aomode;
    pdf,"results/fwhmReconstructionIn" + aomode + "_" + method;
  }else{
    pltitle, meth;
    pdf,"results/statistics/fwhmReconstructionInAllAoModes" + "_" + method;
  }
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

  cnh = readfits("results/seeingProfile.fits");
  l0h = readfits("results/l0hProfile.fits");
  vh = readfits("results/vhProfile.fits");
  alt = readfits("results/altProfile.fits");
  cnh0 = readfits("results/seeing.fits");
  cnhg = readfits("results/groundSeeing.fits");
  cnhalt = readfits("results/altSeeing.fits");
  l0eff = readfits("results/effL0.fits");
  l0g = readfits("results/groundL0.fits");
  l0alt = readfits("results/altL0.fits");
  v0 = readfits("results/effV0.fits");
  v0g = readfits("results/groundV0.fits");
  v0alt = readfits("results/altV0.fits");
  
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
    wa      = where(alt>=20);
    alt(wa) = 20.;
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
func plotScriptErrorBreakdown(Dir,scriptsuffix,path_out,redo=)
/* DOCUMENT eb = plotScriptErrorBreakdown("2013_09_17_onsky/","script293","results/");
 */
{
  local res;
  if(!redo)
    res = readfits(path_out + "scriptProcess_" + scriptsuffix +".fits",err=1);
  
  include, "PRAC_config.i",1;
  dataDirRoot = dataDir + Dir;
  //takes all slopestl files
  pathstl = listVersion(dataDirRoot,"fits","slopestl");
  pathstl = pathstl(where(strfind(scriptsuffix,pathstl)(2,) != -1));
  ntl = numberof(pathstl);
    
  if(is_void(res)){
    s_ir = s_tomo = s_noise = s_bw = s_fit = s_static = s_ncpa = s_ol = [];
    aomode = srsky = srmar = srpar = srborn  = mode = [];
    r0 = L0 = v =  h = ved1 = ved2 = ved3 = ved4 = ved5 = [];
    cnh2 = alt2 = l0h2 = vh2 = []; PSF = array(0.,70,70,ntl);
    res = array(pointer,24);
    
    for(j=1;j<=ntl;j++){
      timetl  = extractDate(pathstl(j));
      h       = grow(h,str2time(strpart(timetl,12:)));
      tmp     = PRAC_main(timetl,Dir = Dir,budgetonly=1);
      mode  = grow(mode,strchar(*tmp(1)));
      
      PSF(,,j)  = *tmp(5);
      
      p       = *tmp(4);
      s_ir  = grow(s_ir,p(1));
      s_fit  = grow(s_fit,p(2));
      s_bw  = grow(s_bw,p(3));
      s_tomo  = grow(s_tomo,p(4));
      s_alias  = grow(s_alias,p(5));
      s_noise  = grow(s_noise,p(6));
      s_ol  = grow(s_ol,p(7));
      s_static = grow(s_static,p(8));
      s_ncpa = grow(s_ncpa,p(9));
      srsky = grow(srsky,p(10));
      srmar = grow(srmar,p(11));
      srpar = grow(srpar,p(12));
      srborn = grow(srborn,p(13));
    
      r0 = grow(r0,(*tmp(2))(1));
      L0 = grow(L0,(*tmp(2))(2));
      v = grow(v,(*tmp(2))(3));

      ved1 = grow(ved1,(*tmp(6))(1));
      ved2 = grow(ved2,(*tmp(6))(2));
      ved3 = grow(ved3,(*tmp(6))(3));
      ved4 = grow(ved4,(*tmp(6))(4));
      ved5 = grow(ved5,(*tmp(6))(5));

      cnh2 = grow(cnh2,(*tmp(3))(,1));
      alt2 = grow(alt2,(*tmp(3))(,2));
      l0h2 = grow(l0h2,(*tmp(3))(,3));
      vh2  = grow(vh2,(*tmp(3))(,4));
      
    }
    res(1) = &mode;res(2) = &PSF;
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

  //Displaying error breakdown
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
  
  labs = ["!s_IR","!s_tomo","!s_alias","!s_n","!s_bw","!s_fit","!s_ol","!s_stat","!s_ncpa"];
  winkill,0;window,0,style="aanda.gs",dpi=90;clr;
  plotsBarDiagram,y1,y2=y2,y3=y3,labs,col1=[char(241)],col2=[char(242)],col3=[char(243)],title=0;
  plt,modes(1),7,.8*max(y1(1),y2(1),y3(1)),tosys=1;
  plt,modes(2),7,.7*max(y1(1),y2(1),y3(1)),tosys=1,color=[128,128,128];
  plt,modes(3),7,.6*max(y1(1),y2(1),y3(1)),tosys=1,color=[64,64,64];
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/" + scriptsuffix +"/" +   scriptsuffix + "_budget.pdf";

  //Displaying averaged PSF

  PSF1 = ((*res(2))(,,1::3))(,,avg);
  PSF2 = ((*res(2))(,,2::3))(,,avg);
  PSF3 = ((*res(2))(,,3::3))(,,avg);
  a=3. ;
  cutmax = (max(max(PSF1),max(PSF2),max(PSF3)))^(1./a);
  box = 70; uz = 0.0297949; 

  window,1;clr;palette,"gray.gp";
  pli,(abs(PSF1))^(1./a),-box*uz/2,-box*uz/2,box*uz/2,box*uz/2,cmin = 0,cmax = cutmax;
  xytitles,"Arcsec","Arcsec";
  pltitle,"On-sky PSF in " + modes(1) + " with SR = " + var2str(arrondi(100*max(PSF1),1)) + "%";
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/" + scriptsuffix + "/"+ scriptsuffix +"_" + modes(1)+ "_psf.pdf";

  window,2;clr;palette,"gray.gp";
  pli,(abs(PSF2))^(1./a),-box*uz/2,-box*uz/2,box*uz/2,box*uz/2,cmin = 0,cmax = cutmax;
  xytitles,"Arcsec","Arcsec";
  pltitle,"On-sky PSF in " + modes(2) + " with SR = " + var2str(arrondi(100*max(PSF2),1)) + "%";
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/" + scriptsuffix + "/" + scriptsuffix +"_" + modes(2)+ "_psf.pdf";

  window,3;clr;palette,"gray.gp";
  pli,(abs(PSF3))^(1./a),-box*uz/2,-box*uz/2,box*uz/2,box*uz/2,cmin = 0,cmax = cutmax;
  xytitles,"Arcsec","Arcsec";
  pltitle,"On-sky PSF in " + modes(3) + " with SR = " + var2str(arrondi(100*max(PSF3),1)) + " %";
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/" + scriptsuffix + "/"+ scriptsuffix +"_" + modes(3)+ "_psf.pdf";

  //Displaying Strehl ratios versus time
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
  
  winkill,4;window,4,style="aanda.gs",dpi=90;clr;
  plmk, SRsky1(sort(h1)),h1(sort(h1)),msize=.4,marker=3,width=10; plg,SRsky1(sort(h1)),h1(sort(h1)),marks=0;
  plmk, SRanal1,h1(sort(h1)),msize=.4,marker=3; plg,SRanal1,h1(sort(h1)),marks=0,type=2;

  plmk, SRsky2(sort(h2)),h2(sort(h2)),msize=.4,marker=4,width=10; plg,SRsky2(sort(h2)),h2(sort(h2)),marks=0;
  plmk, SRanal2,h2(sort(h2)),msize=.4,marker=4; plg,SRanal2,h2(sort(h2)),marks=0,type=2;

  plmk, SRsky3(sort(h3)),h3(sort(h3)),msize=.4,marker=6,width=10; plg,SRsky3(sort(h3)),h3(sort(h3)),marks=0;
  plmk, SRanal3,h3(sort(h3)),msize=.4,marker=6; plg,SRanal3,h3(sort(h3)),marks=0,type=2;
  
  range,0,50;
  xytitles,"Local time (hour)","On-sky SR (%)";
  mm = min(grow(h1,h2,h3))
  plt,"Triangles: " + modes(1),mm,48,tosys=1;
  plt,"Circles: " + modes(2),mm,45,tosys=1;
  plt,"Crosses: " + modes(3),mm,42,tosys=1;
  plt,"Plain line: on-sky SR",mm+.05,48,tosys=1;
  plt,"Dash line: analytic SR",mm+.05,45,tosys=1;
  plt,"Script " + strpart(scriptsuffix,7:),mm+.05,42,tosys=1;
  limits,mm-0.01,mm+0.11;
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/"+ scriptsuffix + "/" + scriptsuffix +"_srvshour.pdf";
  //Displaying VED

   cnh = [median((*res(21))(1::5)),median((*res(21))(2::5)),median((*res(21))(3::5)),median((*res(21))(4::5)),median((*res(21))(5::5))];
   alt = [median((*res(22))(1::5)),median((*res(22))(2::5)),median((*res(22))(3::5)),median((*res(22))(4::5)),median((*res(22))(5::5))]/1000.;
    
  w1 = where(modes!="SCAO");
  if(is_array(w1)){
    m = modes(w1);
    w2 = w1(where(m==m(1)));
    w3 = w1(where(m==m(2)));

    ved1 = sqrt(((*res(20))(w2,)^2)(avg,));
    ved2 = sqrt(((*res(20))(w3,)^2)(avg,));
   
  
    winkill,5;window,5,style="aanda.gs",dpi=90;clr;
    thick=.3;
    for(j=1; j<=5; j++){
      d = 0.5*(0.5-thick);
      d = 2*thick;
      plfp, [char(242)],[0,1,1,0]*ved1(j),([1,1,1,1]*alt(j) -d/2) + [-1,-1,1,1]*thick,[4];
      plfp, [char(244)],[0,1,1,0]*ved2(j),([1,1,1,1]*alt(j) +d/2) + [-1,-1,1,1]*thick,[4];
    }

    xytitles,"Altitude (km)","VED (nm)";
    ml = max(grow(ved1,ved2));
    plt,m(1),1,ml*0.9,tosys=1,color=[128,128,128];
    plt,m(2),1,ml*0.8,tosys=1,color=[64,64,64];
    limits,min(alt)-1,max(max(alt)+1,15);
    range,0,ml*1.05;
    pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/"+ scriptsuffix + "/" + scriptsuffix + "_VED.pdf";
  
  
    //Displaying profiles
   
    
    winkill,6;window,6,style="aanda.gs",dpi=90;clr;
    wmoao = where(modes == "MOAO4L3N");
    suffcnh = readFitsKey(pathstl(wmoao(1)),"CN2H");
    ptrcnh  = restorefits("tomoparam",suffcnh);
    displayLayers,cnh,alt*1000,col=[char(242)],thick=.4,percent=1;
    displayLayers,-(*ptrcnh(2)),*ptrcnh(3),col=[char(241)],percent=1;
    limits,-100,100;
    range,-1,20;
    plt,"Retrieved profile\n on-sky",-90,15,tosys=1;
    plt,"Average profile\n during the script",20,15,tosys=1,color=[128,128,128];
    plt,"Script " + strpart(scriptsuffix,7:),-10,18,tosys=1;
    pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/"+ scriptsuffix + "/" + scriptsuffix + "_profiles.pdf";
  }
  //gives statistics on turbulence
  //r0
  r0 = avg((*res(16))^(-5/3.))^(-3/5.);
  cnh2 = cnh*r0^(-5/3.)/sum(cnh);
  r0g = sum(cnh2(where(alt <= 0.1)))^(-3/5.);
  r0a = sum(cnh2(where(alt >0.1)))^(-3/5.);
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

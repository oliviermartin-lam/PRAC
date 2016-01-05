func define_structs(timedata,verb=)
/* DOCUMENT

 */

{
  extern data;

  restorefits,"slopestl",timedata,pathdata,fake= 1;
  
  data.lambda_vis = 500*1e-9;
  data.nwfs = str2int(readFitsKey(pathdata,"NBWFS"));
  data.its = str2flt(readFitsKey(pathdata,"ITS"));
  data.plscale = str2flt(readFitsKey(pathdata,"PLSCALE"));
  define_tel,pathdata,verb=verb;
  define_wfs,pathdata,verb=verb;
  
  //retrieves engaged data in arcsec
  data.slopes_res = &returnSlopestl(timedata,pathdata,arcsec = 1);
  //retrieves disengaged data in arcsec
  data.slopes_dis = &determinesSlopesdisFromSlopestl(timedata,arcsec=1);
  if(is_scalar(*data.slopes_dis)) return 0;
  data.nframes = dimsof(*data.slopes_dis)(0);
  data.Nslopes = dimsof(*data.slopes_dis)(-1);
  define_dm,pathdata,verb=verb;
  define_rtc,pathdata,verb=verb;
  if(data.wfs(data.its).type == 0){
    return 0;
  }
  define_ir,pathdata,verb=verb;
  define_fourier,verb=verb;
  define_turbu,pathdata,verb=verb;
}

func define_tel(pathdata,verb=)
/* DOCUMENT
 */

{
  extern data;

  data.tel.diam = str2flt(readFitsKey(pathdata,"TELDIAM"));
  data.tel.obs = str2flt(readFitsKey(pathdata,"TELOBS"));
  data.tel.airm = str2flt(readFitsKey(pathdata,"WHTAIRM"));
  data.tel.zen = acos(1./data.tel.airm)*180/pi;
  data.tel.object = strcase(1,readFitsKey(pathdata,"OBJECT"));
}

func define_wfs(pathdata,verb=)
/* DOCUMENT
 */

{

  //WFS
  data.wfs.nssp = readFitsKeyArray(pathdata,"WFSNSLO");
  data.wfs.sX = readFitsKeyArray(pathdata,"WFSSUBX");
  data.wfs.sY = readFitsKeyArray(pathdata,"WFSSUBY");
  data.wfs.type = readFitsKeyArray(pathdata,"WFSTYPE");
  for(i=1;i<=data.nwfs;i++){
    data.wfs(i).pixarc = str2flt(readFitsKey(pathdata,"PIXARC"+int2str(i)));
    data.wfs(i).unit = str2flt(readFitsKey(pathdata,"UNITMRZ"+int2str(i)));
  }
  data.wfs.sym = readFitsKeyArray(pathdata,"WFSSYM");
  //...... Manages LGS ......//
  indlgs = where(data.wfs.type==2);
  indoffngs = readFitsKeyArray(pathdata,"OFFAX");
  nbLgsWfs = numberof(indlgs);
  
  if(is_array(indlgs)){
    tmp = readFitsKey(pathdata,"ALTLGS");
    if(tmp == "ALTLGS not found") tmp = "21000";
    data.wfs(indlgs).lgsH = str2flt(tmp);
    data.wfs(indlgs).lgsRadius = 23;
    data.wfs(indlgs).lgsdH = 1000.;
    data.wfs(indlgs).lgsAngle = 0.0;
      
    // making a square of right diameter : SCIMEASURE CAMERA
    X = [-1,1,-1,1] * data.wfs(indlgs).lgsRadius/sqrt(2);
    Y = [1,1,-1,-1] * data.wfs(indlgs).lgsRadius/sqrt(2);
    // rotating the square
    th = data.wfs(indlgs(1)).lgsAngle * pi/180;
    rotmatrix = [[cos(th),sin(th)],[-sin(th),cos(th)]];
    tmp = rotmatrix(,+)*[X,Y](,+);
    X = tmp(1,);
    Y = tmp(2,);
    
    (data.wfs.x)(indlgs) = X(1:numberof(indlgs));
    (data.wfs.y)(indlgs) = Y(1:numberof(indlgs));
  }

  // recovered positions in mm
  suff_tas = readFitsKey(pathdata,"TAS");
  if(suff_tas == "TAS not found"){
    date = strpart(extractDate(pathdata),strfind("_",extractDate(pathdata))(0)+1:);
    suff_tas = findTruc("tas",date,after=1);
    if(!is_string(suff_tas)){
      suff_tas = findTruc("tas",date);
    }
  }
  tmp = restorefits("tas", suff_tas, pathtas); 
  //tas position the TS is always at the center
  plscale = data.plscale;
  data.wfs(indoffngs).x = tmp(1::2)*plscale;
  data.wfs(indoffngs).y = tmp(2::2)*plscale;
      
  for(i=1;i<=data.nwfs;i++) {
    nssp = data.wfs(i).sX;
    if( nssp==7 ){
      data.wfs(i).mrz = &READFITS("GLOB_mrz_7x7.fits");
      data.wfs(i).zmi = &READFITS("zmi_7x7.fits");
      data.wfs(i).filterouttiltMatrix = &calc_TTMatFilt_1(data.wfs(i).nssp);
      data.wfs(i).keeponlytiltMatrix = &filt_TTvec(data.wfs(i).nssp);
    }
  }
}

func define_dm(pathdata,verb=)
/* DOCUMENT
 */

{
  extern data;

  nssp = data.wfs(data.its).sX;
  
  x = span(-1,1,nssp+1)(,-:1:nssp+1);
  y = transpose(x);
  r = sqrt(x*x+y*y);
  
  msk = r<(1.0+1.4/nssp) & r>(data.tel.obs-1./nssp);
  nn = where( msk );

  data.dm.nactu = numberof(nn);
  data.dm.ndiam = nssp+1;
  data.dm.csX = &(x(nn));
  data.dm.csY = &(y(nn));
  // indices of actuators = [3,4,5,6, 10,11,12,13,14,15, 17 ... ] in a Ndiam X Ndiam image
  data.dm.csI = &(nn);  

  // valeur de x0 a passer a la fonction funcInflu(x,y,x0) pour obtenir un
  // couplage defini pour l'actu voisin, avec x et y exprimes en "rayon pupille"
  couplage = 0.2;
  data.dm.x0 = sqrt(-2/log(couplage)) / nssp;
  
}

func define_rtc(pathdata,verb=)
/* DOCUMENT
 */

{
  extern data;

  //loop
  data.rtc.gain = str2flt(readFitsKey(pathdata,"GAINDM"));
  data.rtc.delay = 0.003;    // loop delay in seconds
  data.rtc.BP = 500.;
  data.rtc.Fe = str2flt(readFitsKey(pathdata,"FREQ"));
  data.rtc.obsmode = strcase(1,readFitsKey(pathdata,"OBS_MODE"));
  data.rtc.rectype = strcase(1,readFitsKey(pathdata,"RECTYPE"));
  wfstype = data.wfs.type
  if(data.rtc.obsmode == "MOAO"){
    suffmt = readFitsKey(pathdata,"MT");
    restorefits,"mt",suffmt,pathmt,fake=1;
    wfstype = readFitsKeyArray(pathmt,"WFSTYPE");
  }
  data.rtc.aomode = giveTomoMode(wfstype,data.rtc.obsmode,data.rtc.rectype);
  data.rtc.ptrlistOffAxisNgs = &readFitsKeyArray(pathdata,"OFFAX");
  data.rtc.ptrlistLgs = &readFitsKeyArray(pathdata,"LGS");
  
  if(data.rtc.obsmode =="MOAO"){
    suffmt = readFitsKey(pathdata,"MT");
    R = restorefits("mt",suffmt);
    //unbiases the reconstructor from sensitivities
    R =  unbiasesMtFromSensitivities(R);
    data.rtc.R = &R;
    if(verb){
      write,format="Calibration time of the reconstructor : %s\n",strpart(suffmt,12:);
    }
  }
  
  
}

func define_ir(pathdata,verb=)
/* DOCUMENT
 */

{
  extern data;
  //Load raw image
  suffir = readFitsKey(pathdata,"IRFILE");
  data.camir.psfraw = &restorefits("ir",suffir,pathir);
  //load bg
  suffbg = readFitsKey(pathir,"BGNAME");
  data.camir.bg = &restorefits("irbg",suffbg);

  if(dimsof(*data.camir.psfraw)(1)==3){
    data.camir.psfraw = data.camir.psfraw(,,avg);
  }
  if(dimsof(data.camir.bg)(1)==3){
    data.camir.bg = data.camir.bg(,,avg);
  }
  //computes the calibrated psf
  data.camir.psfcal = &(*data.camir.psfraw - *data.camir.bg);

  data.camir.camName = strcase(1,readFitsKey(pathir,"IRCAM"));
  tmp = dimsof(*data.camir.psfcal);
  data.camir.dx = tmp(2);
  data.camir.dy = tmp(3);
  data.camir.xPsf = str2int(readFitsKey(pathir,"X0PSF"));
  data.camir.yPsf = str2int(readFitsKey(pathir,"Y0PSF"));
  
  //IR wavelength in m
  data.camir.lambda_ir = str2flt(readFitsKey(pathir,"FILTER"))*1e-9;
  data.camir.uld = str2flt(readFitsKey(pathir,"NPIXLD"));//(D * pixarcIR / lambdaIR)
  data.camir.pixSize = data.camir.uld * data.camir.lambda_ir/data.tel.diam ;
  data.camir.uz = (data.camir.lambda_ir * data.camir.uld/data.tel.diam)*206264.8062471;
  data.camir.deadPixIndex = &CreateDeadPixFrame(readfits("locateDeadPixels"+data.camir.camName+".fits"));
  //size of the cropped on-sky PSF to estimate the Strehl. 
  data.camir.npix_sr = 128;
}

func define_fourier(void,verb= )
/* DOCUMENT

   Allows to create the parameters for doing things in the Fourier
   space, in particular:
   N : size of the support
   uk : size of pixels of the phase spectrum space
   ud : size of pixels of the Dphi space
   dactupix : actuator pitch in ud pixels
   dactu : actuator pitch
   ...
     
   SEE ALSO:
 */
{
  extern data;
  
  nactuX = data.wfs(data.its).sX + 1;
  data.fourier.dactu = data.tel.diam / (nactuX - 1);
  data.fourier.npix = 8 * 2^long(log(nactuX)/log(2)+0.5);
  // number of pixels in the images of phase spectrum (Wfit, Waniso, ...)
  data.fourier.npix = max(data.fourier.npix, 512);

  // one choose the total field of the Dphi image so that its contains a
  // 2*D area with some margin of D/2 at the edges, i.e. a total of 3*D ...

  // taille de l'image de Dphi
  data.fourier.champ_Dphi = 3*data.tel.diam;       

  // size of the pixel of Dphi image (meters)
  data.fourier.ud = data.fourier.champ_Dphi/double(data.fourier.npix);        

  // we now want dactu to be a multiple of ud
  km = long(data.fourier.dactu / data.fourier.ud + 1);
  data.fourier.ud = data.fourier.dactu / km;
  data.fourier.champ_Dphi = data.fourier.npix * data.fourier.ud;
  data.fourier.dactupix = km;

  // size of the pixel of Wiener image (meters^-1)
  data.fourier.uk = 1./data.fourier.champ_Dphi; 
  data.fourier.uz = data.fourier.uk * data.camir.lambda_ir *206264.8062471;
  data.fourier.uld = data.fourier.uk * data.tel.diam;
 
}

func define_turbu(pathdata,verb=)
{
  extern NL_DEF;
  
  data.learn.ttr = 0;

  //determing noise
  data.turbu.varNoise = &determineGlobalParameters(*data.slopes_dis,p,dp,arc=1,verb=verb);
  //global parameters
  data.turbu.r0vis = p(1);
  data.turbu.L0 = p(2);
  data.turbu.v = p(3);
  data.turbu.r0ir = data.turbu.r0vis * (data.camir.lambda_ir/data.lambda_vis)^1.2;
  //Concatenating uncertainties on global parameters into data.uncertainties struct
  data.uncertainties.r0  = dp(1);
  data.uncertainties.L0  = dp(2);
  data.uncertainties.v   = dp(3);

  if(verb){
    write,format="r0(%g nm)     = %.3g cm +/- %.2g\n",data.lambda_vis*1e9,data.turbu.r0vis*100., 100*dp(1);
    write,format="L0             = %.3g m +/- %.2g\n",data.turbu.L0,dp(2);
    write,format="Windspeed      = %.3g m/s +/- %.2g\n",data.turbu.v,dp(3);
  }
  
  //loading of the on-sky retrieved profiles from post-processing
  data.learn.nl = NL_DEF;
  
  goodDir = "profiles/" + procDir;
  ptr_prof = readfits(goodDir + "profiles_" + extractDate(pathdata) + "_nl_" + var2str(NL_DEF) +".fits",err=1,verb=verb);

  if(is_pointer(ptr_prof) && verb) write,"\rReading of the profiles successfully done";
    
  if(!is_pointer(ptr_prof)){
    //getting the cnh and l0h profiles
    ptrcn2h = fitCovarianceMatrix(*data.slopes_dis,NL_DEF,fitl0=2,FitMethod=2,fullHD=0,tomores=20000./(NL_DEF+1),ttr=1,verb=verb);
    //getting the windspeed profiles
    getWindspeedProfile,dvh,ddirh,verb=verb;
    //merging all retrieved parameters and uncertainties
    ptr_prof = array(pointer,12);
    for(i=1;i<=4;i++){ptr_prof(i) = ptrcn2h(i);}
    ptr_prof(5)  = &data.learn.vh;
    ptr_prof(6)  = &data.learn.dirh;
    ptr_prof(7)  = &data.uncertainties.cnh;
    ptr_prof(8)  = &data.uncertainties.altitude;
    ptr_prof(9)  = &data.uncertainties.l0h;
    ptr_prof(10) = &data.uncertainties.tracking;
    ptr_prof(11) = &data.uncertainties.vh;
    ptr_prof(12) = &data.uncertainties.dirh;
    
    //writing the files
    if(!direxist(goodDir)) system,"mkdir " + goodDir;
    writefits, goodDir + "profiles_" + extractDate(pathdata) + "_nl_" + var2str(NL_DEF)+".fits",ptr_prof;
  }

  //cnh profile
  data.learn.cnh               = abs(*ptr_prof(1));
  data.uncertainties.cnh       = abs(*ptr_prof(7));
  data.learn.altitude          = *ptr_prof(2);
  data.uncertainties.altitude  = abs(*ptr_prof(8));
  data.learn.l0h               = abs(*ptr_prof(3));
  data.uncertainties.l0h       = abs(*ptr_prof(9));
  data.learn.tracking          = *ptr_prof(4);
  data.uncertainties.tracking  = abs(*ptr_prof(10));

  getWindspeedProfile,dvh,ddirh,verb=verb;
  ptr_prof(5)  = &data.learn.vh;
  ptr_prof(11) = &data.uncertainties.vh;
  
  data.learn.vh                = abs(*ptr_prof(5));
  data.uncertainties.vh        = abs(*ptr_prof(11));


  w = where(data.learn.l0h >=100.);
 
}
  

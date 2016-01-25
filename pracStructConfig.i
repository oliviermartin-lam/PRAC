func define_structs(timedata,simu=,verb=)
/* DOCUMENT define_structs,"02h30m48s",verb=1;

 */

{
  extern ast,tel,atm,rtc,wfs,dm,cam,sys,learn,learn_fit,covMatrix,budget,otf,psf;
  if(simu == 1){
    pracRead,pathdata;
    return 0;
  }else{
    restorefits,"slopestl",timedata,pathdata,fake = 1;
    nwfs = str2int(readFitsKey(pathdata,"NBWFS"));
    nngs = str2int(readFitsKey(pathdata,"NOFFAX")) + 1;
    nl   = 5;
  
    //Structures instantiation
    include, "pracStructInit.i",1;
    ast    = asterism_struct();
    tel    = telescope_struct();
    atm    = atm_struct();
    rtc    = rtc_struct();
    wfs    = array(wfs_struct,nwfs);
    dm     = dm_struct();
    cam    = cam_struct();
    sys    = sys_struct();
    learn  = learn_struct();
    learn_fit  = learnfit_struct();
    covMatrix = covmatrix_struct();
    budget = budget_struct();
    otf    = otf_struct();
    psf    = psf_struct();

    // Structures filling
    defineAsterism,pathdata,verb=verb;
    defineCam,pathdata,verb=verb;
    defineTel,pathdata,verb=verb;
    defineWfs,pathdata,verb=verb;
    defineRtc,timedata,verb=verb;
    defineDm,pathdata,verb=verb;
    defineSys,tracking;
    defineAtm,pathdata,varNoise,verb=verb;
    defineLearn;
    defineCovMat,varNoise,verb=verb;
  }
}

func defineAsterism(pathdata,verb=)
/* DOCUMENT
 */
{
  extern ast;
  ast.name = strcase(1,readFitsKey(pathdata,"OBJECT"));

  if(ast.name=="A47" | ast.name=="AST47"){
    alpha  = "21 12 3.0";
    delta = "38 36 45.0";
    magAstInR = [8.7,10.2,9.9,11.]
  }
  if(ast.name=="A453" | ast.name=="AST453"){
    alpha  = "+1 36 44.400";
    delta = "47 21 23.0";
    magAstInR = [11.79,10.76,9.92,10.92];
  }
  if(ast.name=="A396" | ast.name=="AST396"){
    alpha  = "6 48 15.0";
    delta = "41 05 44.0";
    magAstInR = [8.988,10.875,8.902];
  }
  if(ast.name=="A53" | ast.name=="AST53"){
    alpha  = "+23 24 30.0";
    delta = "40 53 55.0";
    magAstInR = [10.9,11.2,9.9,9.8];
  }
 if(ast.name=="A110" | ast.name=="AST110"){
    alpha  = "+19 31 11";
    delta = "34 56 02";
    magAstInR = [10.9,9.7,11.0,11.6];
  }
 if(ast.name=="AS2" | ast.name=="A32" | ast.name=="AST32"){
    alpha  = "+18 27 7.0";
    delta = "26 52 34.0";
    magAstInR = [12.1,11.2,11.7,11.75];
  }
 if(ast.name=="A34" | ast.name=="AST34"){
    alpha  = "+18 51 30.0";
    delta = "10 20 3.0";
    magAstInR = [9.7,9.2,9.7,11.6];
  }
 if(ast.name=="A51" | ast.name=="AST51"){
    alpha  = "+22 48 5";
    delta = "39 17 3.0";
    magAstInR = [11.3,8.2,11.6,9.9];
  }

 if(ast.name=="A10" | ast.name=="AST10"){
    alpha  = "+5 52 17.0";
    delta = "32 34 36.0";
    magAstInR = [11.4,9.4,11.6,10.8];
  }
 
 if(ast.name=="AT1" | ast.name=="ASTT1"){
    alpha  = "+13 41 38.2";
    delta = "+7 36 21.0";
    magAstInR = [9.0,0,9.8,9.6];
  }

 if(ast.name=="A12" | ast.name=="AST12"){
    alpha  = "+6 1 9.0";
    delta = "+23 20 29.0";
    magAstInR = [8.3,11.2,10.7,10.0];
  }

 
 ast.nOffAxisSources = str2flt(readFitsKey(pathdata,"NOFFAX"));
 ast.rightAscension = alpha;
 ast.declinaison = delta;
 ast.magnitude.R = magAstInR;
 ast.photometry  = str2flt(readFitsKey(pathdata,"FILTER"))*1e-9;

}

func defineCam(pathdata,verb=)
/* DOCUMENT
 */

{
  extern cam;

  //size of the cropped on-sky PSF to estimate the Strehl. 
  cam.nPixelsCropped = 128;
    
  //Load raw image
  suffir = readFitsKey(pathdata,"IRFILE");
  if(suffir != "IRFILE not found"){
    cam.psfRaw = &restorefits("ir",suffir,pathir);
    //load bg
    suffbg = readFitsKey(pathir,"BGNAME");
    cam.bg = &restorefits("irbg",suffbg);

    if(dimsof(*cam.psfRaw)(1)==3){
      cam.psfRaw = cam.psfRaw(,,avg);
    }
    if(dimsof(cam.bg)(1)==3){
      cam.bg = cam.bg(,,avg);
    }

    cam.name = strcase(1,readFitsKey(pathir,"IRCAM"));
    tmp = dimsof(*cam.psfRaw);
    cam.nxPixels = tmp(2);
    cam.nyPixels = tmp(3);
    cam.xPsf = str2int(readFitsKey(pathir,"X0PSF"));
    cam.yPsf = str2int(readFitsKey(pathir,"Y0PSF"));
    cam.exposureTime = str2flt(readFitsKey(pathdata,"EXPTIME"))*1e-6;//in sec

    //IR wavelength in m
    cam.lambda = str2flt(readFitsKey(pathir,"FILTER"))*1e-9;
    cam.dPixSizeOnLambda = str2flt(readFitsKey(pathir,"NPIXLD"));
    D = str2flt(readFitsKey(pathdata,"TELDIAM"));
    cam.pixSize = cam.dPixSizeOnLambda * cam.lambda/D*radian2arcsec ;
    cam.deadPixIndex = &CreateDeadPixFrame(readfits("fitsFiles/locateDeadPixels"+cam.name+".fits"));
    //computes the calibrated psf

    cam.psfCal = &processCanaryPSF(suffir,box=cam.nPixelsCropped);
    otf.sky    = &(roll(fft(roll(*cam.psfCal)).re));
  }
}

func defineTel(pathdata,verb=)
/* DOCUMENT
 */

{
  extern tel;
  tel.diam   = str2flt(readFitsKey(pathdata,"TELDIAM"));
  tel.obs    = str2flt(readFitsKey(pathdata,"TELOBS"));
  airm       = readFitsKey(pathdata,"WHTAIRM");
  if(airm != "WHTAIRM not found"){
    tel.airmass   = str2flt(airm);
    tel.zenith    = acos(1./tel.airmass)*180/pi;
  }
  
  tel.plateScale = str2flt(readFitsKey(pathdata,"PLSCALE"));
  tel.azimuth = 0;

  // number of pixels in the images of phase spectrum (Wfit, Waniso, ...)
  nLenslet = max(readFitsKeyArray(pathdata,"WFSSUBX"));
  tel.pitch =  tel.diam / nLenslet;
  tel.nPixels = max(8*2^long(log(nLenslet+1)/log(2)+0.5),512);
  tel.pixSize = cam.pixSize
  // field of View in meters
  tel.foV = tel.nPixels*tel.pixSize;
  //actuator pitch in nPixels in forcing dm.dactu as a tel.pixSize multiple
  tel.pitchSize = tel.pitch / tel.pixSize; 
  // pixel size in the Fourier domain
  tel.fourierPixSize = 1./tel.foV;
  tel.lambdaPerFov = tel.fourierPixSize * ast.photometry * radian2arcsec; 
  tel.dPerFov = tel.fourierPixSize * tel.diam;
  tel.airyPeak = (pi/4)*(1 - tel.obs^2) * cam.dPixSizeOnLambda^2;

}

func defineWfs(pathdata,verb=)
/* DOCUMENT
 */

{
  extern wfs;
  //WFS
  wfs.nValidMeas(1:nwfs)  = readFitsKeyArray(pathdata,"WFSNSLO");
  wfs.nValidSubAp(1:nwfs) = wfs.nValidMeas/2;
  wfs.nLenslet(1:nwfs)    = readFitsKeyArray(pathdata,"WFSSUBX");
  wfs.type(1:nwfs)        = readFitsKeyArray(pathdata,"WFSTYPE");
  wfs.nPixels(1:nwfs)     = 16;
  
  for(i=1;i<=nwfs;i++){
    wfs(i).pixSize = str2flt(readFitsKey(pathdata,"PIXARC"+int2str(i)));
    wfs(i).unit    = str2flt(readFitsKey(pathdata,"UNITMRZ"+int2str(i)));
  }
  wfs.sym(1:nwfs) = readFitsKeyArray(pathdata,"WFSSYM");


  // Recovering off-axis stars position from TAS calibration
  indOffNgs = readFitsKeyArray(pathdata,"OFFAX");
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
  wfs(indOffNgs).x = tmp(1::2)*tel.plateScale;
  wfs(indOffNgs).y = tmp(2::2)*tel.plateScale;

  
  //...... Managing LGS ......//
  indLgs = where(wfs.type==2);
  nLgsWfs = numberof(indLgs);

  if(is_array(indLgs)){
    if(numberof(indLgs) == 1){
      wfs(indLgs).lgsH = 13500.;
    }else{
      tmp = readFitsKey(pathdata,"ALTLGS");
      if(tmp == "ALTLGS not found") tmp = "21000";
      wfs(indLgs).lgsH = str2flt(tmp);
      wfs(indLgs).lgsRadius = 23;
      wfs(indLgs).lgsdH = 1000.;
      wfs(indLgs).lgsAngle = 0.0;
      
      // making a square of right diameter : SCIMEASURE CAMERA
      X = [-1,1,-1,1] * wfs(indLgs).lgsRadius/sqrt(2);
      Y = [1,1,-1,-1] * wfs(indLgs).lgsRadius/sqrt(2);
      // rotating the square
      th = wfs(indLgs(1)).lgsAngle * pi/180;
      rotMatrix = [[cos(th),sin(th)],[-sin(th),cos(th)]];
      tmp = rotMatrix(,+)*[X,Y](,+);
      X = tmp(1,);
      Y = tmp(2,);
    
      wfs(indLgs).x = X(1:numberof(indLgs));
      wfs(indLgs).y = Y(1:numberof(indLgs));
    }
  }
}

func defineRtc(timedata,verb=)
/* DOCUMENT
 */

{
  extern rtc;

  // .... System parameters retrieval ....
  rtc.nWfs = nwfs;
  rtc.its     = str2int(readFitsKey(pathdata,"ITS")); 
  rtc.loopGain = str2flt(readFitsKey(pathdata,"GAINDM"));
  rtc.BP = 500.;
  rtc.Fe = str2flt(readFitsKey(pathdata,"FREQ"));
  rtc.delay = 0.003;    // loop delay in seconds
  rtc.frameDelay = rtc.delay*rtc.Fe + 1.05;
  rtc.obsMode = strcase(1,readFitsKey(pathdata,"OBS_MODE"));
  rtc.nWfs = nwfs;
  rtc.ptrListOffAxisNgs = &readFitsKeyArray(pathdata,"OFFAX");
  rtc.ptrListLgs = &readFitsKeyArray(pathdata,"LGS");
  rtc.nLgsWfs = numberof(rtc.ptrListLgs);
  rtc.nNgsWfs = numberof(rtc.ptrListOffAxisNgs)+1;


  // .... matrices retrieval ....
  //interation matrix
  suffmi  = readFitsKey(pathdata,"MI");//in pixel/ADU
  mi = volt2adu*wfs(rtc.its).pixSize*restorefits("mi",suffmi);//in arcsec/V
  rtc.mi = &(mi(,1:-2)); 
  //command matrix
  suffmc  = readFitsKey(pathdata,"MC");//in ADU/pixel
  mc = restorefits("mc",suffmc)/(volt2adu*wfs(rtc.its).pixSize);
  rtc.mc = &(mc(1:-2,));
  
  if(is_void(mc)){
    rtc.mi = &intermat(dm);
    rtc.mc = &computeCommandMat(mi, nmf=-1, condi=30., disp=0);
  }
  
  // .... RTC provided vectors
  //voltages
  suffvolts = readFitsKey(pathdata,"VOLTFILE");
  if(suffvolts != "VOLTFILE not found"){
    ptr_volts = restorefits("voltstl",suffvolts,pathvolts);
    if(dimsof(*ptr_volts(1))(1) != 0)
      rtc.volts =&((adu2volt*(*ptr_volts(1)))(1:-2,));//exclusion of steering mirror for LGS
  }

  // .... Residual slopes
  rtc.slopes_res = &returnSlopestl(timedata,pathdata,arcsec = 1);

  // .... Reconstructed uncorrected slopes
  rtc.slopes_dis = &determinesSlopesdisFromSlopestl(timedata,arcsec=1);
  rtc.nFrames = dimsof(*rtc.slopes_dis)(0);
  rtc.nSlopes = dimsof(*rtc.slopes_dis)(-1);

  // ..... MOAO Reconstructor
  if(rtc.obsMode == "MOAO"){
    suffmt = readFitsKey(pathdata,"MT");
    restorefits,"mt",suffmt,pathmt,fake=1;
    wfs.type = readFitsKeyArray(pathmt,"WFSTYPE");
    R = restorefits("mt",suffmt);
    //unbiases the reconstructor from sensitivities
    R =  unbiasesMtFromSensitivities(R);
    rtc.R = &R;
    if(verb){
      write,format="Calibration time of the reconstructor : %s\n",strpart(suffmt,12:);
    }
  }
  rtc.recType = strcase(1,readFitsKey(pathdata,"RECTYPE"));
  rtc.aoMode = giveTomoMode(wfs.type,rtc.obsMode,rtc.recType);
}

func defineAtm(pathdata,&varNoise,verb=)
{
  extern atm;
  atm.lambda = 500e-9;
  //determing noise
  varNoise = determineGlobalParameters(*rtc.slopes_dis,p,dp,arc=1,verb=verb);
  //global parameters
  atm.r0 = p(1);
  atm.L0 = p(2);
  atm.v =  p(3);
  
  if(verb){
    write,format="r0(%g nm)     = %.3g cm +/- %.2g\n",atm.lambda*1e9,atm.r0*100., 100*dp(1);
    write,format="L0             = %.3g m +/- %.2g\n",atm.L0,dp(2);
    write,format="Windspeed      = %.3g m/s +/- %.2g\n",atm.v,dp(3);
  }


  //noise covariance matrix
  covnoise = array(0.,rtc.nSlopes,rtc.nSlopes);
  takesDiag(covnoise) = varNoise;
  covnoise = handle_tilt_from_wfstype(covnoise);
  covMatrix.noise = &covnoise;
  
  //loading of the on-sky retrieved profiles from post-processing
  atm.nLayers = nl;

  if(strpart(procDir,0:0) != "/")
    procDir = procDir+"/";
  
  goodDir = "profiles/" + procDir;
  ptr_prof = readfits(goodDir + "profiles_" + extractDate(pathdata) + "_nl_" + var2str(nl) +".fits",err=1,verb=verb);

  if(is_pointer(ptr_prof) && verb) write,"\rReading of the profiles successfully done";
    
  if(!is_pointer(ptr_prof)){
    //getting the cnh and l0h profiles
    ptrcn2h = fitCovarianceMatrix(*rtc.slopes_dis,nl,fitl0=2,FitMethod=2,fullHD=0,tomores=20000./(nl+1),ttr=1,verb=verb);
    //getting the windspeed profiles
    getWindspeedProfile,dvh,ddirh,verb=verb;
    //merging all retrieved parameters and uncertainties
    ptr_prof = array(pointer,6);
    for(i=1;i<=4;i++){ptr_prof(i) = ptrcn2h(i);}
    ptr_prof(5)  = &atm.vh;
    ptr_prof(6)  = &atm.dirh;
    
    //writing the files
    if(!direxist(goodDir)) system,"mkdir " + goodDir;
    writefits, goodDir + "profiles_" + extractDate(pathdata) + "_nl_" + var2str(nl)+".fits",ptr_prof;
  }

  //Updating profiles
  atm.cnh      = abs(*ptr_prof(1));
  atm.altitude = *ptr_prof(2);
  atm.l0h      = abs(*ptr_prof(3));
  atm.vh       = abs(*ptr_prof(5));
  sys.tracking     = *ptr_prof(4);
}


func defineDm(pathdata,verb=)
/* DOCUMENT
 */

{
  extern dm;
  nLenslet = wfs(rtc.its).nLenslet;
  dm.pitch = tel.diam / nLenslet;
  dm.nActu = nLenslet+1;
  //retrieving the number of valid actuators in the pupil
  x = span(-1,1,nLenslet+1)(,-:1:nLenslet+1);
  y = transpose(x);
  r = sqrt(x*x+y*y);
  msk = r<(1.0+1.4/nLenslet) & r>(tel.obs-1./nLenslet);
  nn = where( msk );
  dm.nValidActu = numberof(nn);
  dm.csX = &(x(nn));
  dm.csY = &(y(nn));
  // indices of actuators = [3,4,5,6, 10,11,12,13,14,15, 17 ... ] in a Ndiam X Ndiam image
  dm.csI = &(nn);  

  // valeur de x0 a passer a la fonction funcInflu(x,y,x0) pour obtenir un
  // couplage defini pour l'actu voisin, avec x et y exprimes en "rayon pupille"
  dm.couplage = 0.2;
  dm.x0 = sqrt(-2/log(dm.couplage)) / nLenslet;
  
}


func defineSys(void)
/* DOCUMENT
 */
{
  extern sys;
  sys.xshift = array(0.,nwfs);
  sys.yshift = array(0.,nwfs);
  sys.theta  = array(0.,nwfs); 
  sys.magnification = array(1.,nwfs);
  sys.centroidGain = array(1.,nwfs);
  sys.slopesToZernikeMatrix = &readfits("fitsFiles/GLOB_mrz_7x7.fits");
  sys.zernikeToSlopesMatrix = &readfits("fitsFiles/zmi_7x7.fits");
}
 
func defineLearn(void)
/*DOCUMENT
 */
{
  extern learn;

  learn.nl = atm.nLayers;
  learn.ttr = 0;
  learn.diagonal = 0;
  learn.runLearnTwoSteps = 0;
  learn.cnh = atm.cnh;
  learn.altitude = atm.altitude;
  learn.l0h = atm.l0h;
  learn.tracking = sys.tracking;
  learn.xshift = sys.xshift;
  learn.yshift = sys.yshift;
  learn.magnification = sys.magnification;
  learn.theta = sys.theta;
  learn.centroidGain = sys.centroidGain;
}

func defineCovMat(varNoise,verb=)
/*DOCUMENT
 */
{
  extern covMatrix;

  //slopes covariance matrix
  if(verb){
      write,"\rComputing slopes covariance matrix...";
  }
  s = *rtc.slopes_dis;
  s -= s(,avg);
  covslopes = s(,+)*s(,+)/rtc.nFrames;
  covslopes  = handle_tilt_from_wfstype(covslopes);
  covMatrix.slopes = &covslopes;
  

  //Phase gradient synthetic matrix
  if(verb){
      write,"\rComputing synthetic slopes covariance matrix...";
  }
  fitEstim = packcoeffs(learn);
  covMatrix.learn = &covMatModel(learn, fitEstim,loworder=0,verb=verb);

  //Contribution of correctable modes by the AO system into cthe Phase gradient matrix
  if(verb){
      write,"\rComputing parallel modes covariance matrix...";
  }
  covMatrix.parallel = &covMatModel(learn, fitEstim,loworder=1,verb=verb);

  //Contribution of correctable modes by the Ao system into
  if(verb){
      write,"\rComputing aliasing covariance matrix...";
  }
  covMatrix.aliasing = &(*covMatrix.learn - *covMatrix.parallel);

  // Isoplanatic contribution from something else turbulence (vibration, tracking telescope...)
  covTracking = 0*(*covMatrix.learn);
  covMatrix.tracking = &trackingMatCov(sys.tracking, covTracking);
}

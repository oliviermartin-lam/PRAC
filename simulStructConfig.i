func defineSimulStructs(void,verb=)
/* DOCUMENT defineSimulStructs

 */

{
  extern ast,tel,atm,rtc,wfs,sys,cam,dm,learn,learn_fit,covMatrix;

  nwfs = 8;nngs=4;nl=4;

  include, "pracStructInit.i",1;
  ast    = asterism_struct();
  cam    = cam_struct();
  tel    = telescope_struct();
  atm    = atm_struct();
  rtc    = rtc_struct();
  wfs    = array(wfs_struct,nwfs);
  dm     = dm_struct();    
  sys    = sys_struct();
  cam    = cam_struct();
  dm     = dm_struct();
  learn  = learn_struct();
  learn_fit = learnfit_struct();
  covMatrix = covmatrix_struct();
  
  /*
     ____   ___  _   _ ____   ____ _____ ____  
    / ___| / _ \| | | |  _ \ / ___| ____/ ___| 
    \___ \| | | | | | | |_) | |   |  _| \___ \ 
     ___) | |_| | |_| |  _ <| |___| |___ ___) |
    |____/ \___/ \___/|_| \_\\____|_____|____/ 
                                           
  */

  ast.nOffAxisSources = 3;
  //ast.magnitude.R = [10,10,10,10];
  ast.photometry  = 1.66e-6;
  wfs.x = [-16.2635,16.2635,-16.2635,16.2635,-43.4014,-33.3075,46.3348,0];
  wfs.y = [16.2635,16.2635,-16.2635,-16.2635,-29.2985,24.7061,-11.6072,0];

  /*
   _____ _____ _     _____ ____   ____ ___  ____  _____ 
  |_   _| ____| |   | ____/ ___| / ___/ _ \|  _ \| ____|
    | | |  _| | |   |  _| \___ \| |  | | | | |_) |  _|  
    | | | |___| |___| |___ ___) | |__| |_| |  __/| |___ 
    |_| |_____|_____|_____|____/ \____\___/|_|   |_____|
                                                      
  */

  tel.diam    = 4.2;
  tel.obs     = 0.285;
  tel.airmass = 1.;
  tel.zenith  = acos(1./tel.airmass)*180/pi;  
  tel.azimuth = 0;

  // number of pixels in the images of phase spectrum (Wfit, Waniso, ...)
  nLenslet = 32;
  tel.pitch =  tel.diam / nLenslet;
  tel.nPixels = max(8*2^long(log(nLenslet+1)/log(2)+0.5),512);
  tel.pixSize = 0.03;
  // field of View in meters
  tel.foV = tel.nPixels*tel.pixSize;
  //actuator pitch in nPixels in forcing dm.dactu as a tel.pixSize multiple
  tel.pitchSize = tel.pitch / tel.pixSize; 
  // pixel size in the Fourier domain
  tel.fourierPixSize = 1./tel.foV;
  tel.lambdaPerFov = tel.fourierPixSize * ast.photometry * radian2arcsec; 
  tel.dPerFov  = tel.fourierPixSize * tel.diam;
  tel.aera     = tel.diam^2*(1-tel.obs^2)*pi/4;
  tel.aeraInPix= tel.aera / tel.pixSize^2;

  /*
 __        _______ ____  
 \ \      / /  ___/ ___| 
  \ \ /\ / /| |_  \___ \ 
   \ V  V / |  _|  ___) |
    \_/\_/  |_|   |____/ 

  */
  
  wfs.nValidMeas(1:nwfs)  = 72;
  wfs.nValidSubAp(1:nwfs) = wfs.nValidMeas/2;
  wfs.nLenslet(1:nwfs)    = 7;
  wfs.type(1:nwfs)        = [2,2,2,2,1,1,1,1];
  wfs.nPixels(1:nwfs)     = 16;
  wfs(i).pixSize	  = 0.025;
  wfs.sym(1:nwfs)         = 1;

  
  //...... Managing LGS ......//
  indLgs = where(wfs.type==2);
  nLgsWfs = numberof(indLgs);

  if(is_array(indLgs)){
    wfs(indLgs).lgsH = 21000;
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
  

  /*
     ____ _____ ____ 
    |  _ \_   _/ ___|
    | |_) || || |    
    |  _ < | || |___ 
    |_| \_\|_| \____|
                 
  */
  rtc.nWfs     = nwfs;
  rtc.its      = rtc.nWfs;
  rtc.loopGain = 1.;
  rtc.BP       = 10e3;
  rtc.Fe       = 150.; 
  rtc.delay    = 0.003;    // loop delay in seconds
  rtc.frameDelay = rtc.delay*rtc.Fe + 1.05;
  rtc.obsMode = "MOAO";

  rtc.ptrListOffAxisNgs = &[5,6,7];
  rtc.ptrListLgs = &[1,2,3,4];
  rtc.nLgsWfs    = numberof(*rtc.ptrListLgs);
  rtc.nNgsWfs    = numberof(*rtc.ptrListOffAxisNgs)+1;
  rtc.nSlopes    = sum(wfs.nValidMeas);
  
  /*
        _  _____ __  __  ___  ____  ____  _   _ _____ ____  _____ 
       / \|_   _|  \/  |/ _ \/ ___||  _ \| | | | ____|  _ \| ____|
      / _ \ | | | |\/| | | | \___ \| |_) | |_| |  _| | |_) |  _|  
     / ___ \| | | |  | | |_| |___) |  __/|  _  | |___|  _ <| |___ 
    /_/   \_\_| |_|  |_|\___/|____/|_|   |_| |_|_____|_| \_\_____|

  */
  atm.r0 = 0.15;	
  atm.L0 = 10.;
  atm.v =  10.;
  atm.nLayers  = nl;
  atm.cnh      = atm.r0^(-5/3.) * [65,10,15,10]/100.;
  atm.altitude = [0.,4.,10.,15.]*1e3;
  atm.l0h      = [10.,10.,10.,10.];
  atm.vh       = [8.,12.,15.,18.];
  atm.lambda   = 500e-9;

    
  /*
     ______   ______ _____ _____ __  __ 
    / ___\ \ / / ___|_   _| ____|  \/  |
    \___ \\ V /\___ \ | | |  _| | |\/| |
     ___) || |  ___) || | | |___| |  | |
    |____/ |_| |____/ |_| |_____|_|  |_|
                                    
  */

  sys.tracking = [0.,0.,0.];
  sys.xshift = array(0.,nwfs);
  sys.yshift = array(0.,nwfs);
  sys.theta  = array(0.,nwfs); 
  sys.magnification = array(1.,nwfs);
  sys.centroidGain = array(1.,nwfs);

  

  /*
     ____  __  __ 
    |  _ \|  \/  |
    | | | | |\/| |
    | |_| | |  | |
    |____/|_|  |_|
              
  */
  
  dm.pitch = tel.diam / nLenslet;
  dm.nActu = nLenslet+1;
  //retrieving the number of valid actuators in the pupil in a Fried's geometry
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
  
/*
 ____   ____ ___ _____ _   _  ____ _____ 
/ ___| / ___|_ _| ____| \ | |/ ___| ____|
\___ \| |    | ||  _| |  \| | |   |  _|  
 ___) | |___ | || |___| |\  | |___| |___ 
|____/ \____|___|_____|_| \_|\____|_____|
                                         
  ____    _    __  __ _____ ____      _    
 / ___|  / \  |  \/  | ____|  _ \    / \   
| |     / _ \ | |\/| |  _| | |_) |  / _ \  
| |___ / ___ \| |  | | |___|  _ <  / ___ \ 
 \____/_/   \_\_|  |_|_____|_| \_\/_/   \_\
                                           
*/

  cam.lambda = ast.photometry;

  /*
 _     _____    _    ____  _   _ 
| |   | ____|  / \  |  _ \| \ | |
| |   |  _|   / _ \ | |_) |  \| |
| |___| |___ / ___ \|  _ <| |\  |
|_____|_____/_/   \_\_| \_\_| \_|
                                 
  */
  
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

  
  /*
 ____  _____ ____ ___  _   _ ____ _____ ____  _   _  ____ _____ ___  ____  
|  _ \| ____/ ___/ _ \| \ | / ___|_   _|  _ \| | | |/ ___|_   _/ _ \|  _ \ 
| |_) |  _|| |  | | | |  \| \___ \ | | | |_) | | | | |     | || | | | |_) |
|  _ <| |__| |__| |_| | |\  |___) || | |  _ <| |_| | |___  | || |_| |  _ < 
|_| \_\_____\____\___/|_| \_|____/ |_| |_| \_\\___/ \____| |_| \___/|_| \_\
                                                                           
  */

  //Atmosphere covariance matrix
 
  fitEstim = packcoeffs(learn);
  covMatrix.atm = &covMatModel(learn,fitEstim,verb=verb);

  // Noise covariance matrix
  covnoise = array(0.,rtc.nSlopes,rtc.nSlopes);
  covMatrix.noise = &covnoise;

  
  //Introducing mis-registration
  learn.nl = nl;
  //learn.cnh(2:) = atm.cnh(2:)*0.;
  //learn.l0h(1:) = 100;
  //learn.xshift(1:-1) = .1;
  //learn.altitude = [0,2000,8000,18000];
  fitEstim = packcoeffs(learn);
  
  covMatrix.learn    = &covMatModel(learn,fitEstim,verb=verb);
  covMatrix.parallel = &covMatModel(learn,fitEstim,loworder=1,verb=verb);
  covMatrix.aliasing = &(*covMatrix.learn - *covMatrix.parallel);
  covTracking        = 0*(*covMatrix.learn);
  covMatrix.tracking = &trackingMatCov(learn.tracking, covTracking);

  //Interaction matrix
  rtc.mi = &intermat(dm);
  rtc.mc = &computeCommandMat(*rtc.mi, nmf=-1, condi=30., disp=0);

  //Tomographic reconstructor
  Coffoff = (*covMatrix.learn +  *covMatrix.noise)(norange(rtc.its),norange(rtc.its));
  Conoff  = (*covMatrix.learn)(slrange(rtc.its),norange(rtc.its));

  //inversion
  iCoffoff = invgen(Coffoff,1,cond=100);
  rtc.R    = &(Conoff(,+) * iCoffoff(+,));
}






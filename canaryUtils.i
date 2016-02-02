/*
    _    ____ _____ _____ ____  ___ ____  __  __ 
   / \  / ___|_   _| ____|  _ \|_ _/ ___||  \/  |
  / _ \ \___ \ | | |  _| | |_) || |\___ \| |\/| |
 / ___ \ ___) || | | |___|  _ < | | ___) | |  | |
/_/   \_\____/ |_| |_____|_| \_\___|____/|_|  |_|
                                                 
*/

func givesCoordinatesFromObject(object)
{
  object = strcase(1,object);

  if(object=="A47" | object=="AST47"){
    alpha  = "21 12 3.0";
    delta = "38 36 45.0";
  }
  if(object=="A453" | object=="AST453"){
    alpha  = "+1 36 44.400";
    delta = "47 21 23.0";
  }
  if(object=="A396" | object=="AST396"){
    alpha  = "6 48 15.0";
    delta = "41 05 44.0";
  }
  if(object=="A53" | object=="AST53"){
    alpha  = "+23 24 30.0";
    delta = "40 53 55.0";
  }
 if(object=="A110" | object=="AST110"){
    alpha  = "+19 31 11";
    delta = "34 56 02";
  }
 if(object=="AS2" | object=="A32" | object=="AST32"){
    alpha  = "+18 27 7.0";
    delta = "26 52 34.0";
  }
 if(object=="A34" | object=="AST34"){
    alpha  = "+18 51 30.0";
    delta = "10 20 3.0";
  }
 if(object=="A51" | object=="AST51"){
    alpha  = "+22 48 5";
    delta = "39 17 3.0";
  }

 if(object=="A10" | object=="AST10"){
    alpha  = "+5 52 17.0";
    delta = "32 34 36.0";
  }
 
 if(object=="AT1" | object=="ASTT1"){
    alpha  = "+13 41 38.2";
    delta = "+7 36 21.0";
  }

 if(object=="A12" | object=="AST12"){
    alpha  = "+6 1 9.0";
    delta = "+23 20 29.0";
  }
 
 
  return [alpha,delta];

}

func airmassFromDate(date,alphaJ2000,deltaJ2000)
/* DOCUMENT
 */
{
  //extract date from path
  Y = str2int(strpart(date,1:4));
  M = str2int(strpart(date,6:7));
  D = str2int(strpart(date,9:10));
  H = str2int(strpart(date,12:13))-1.;//correction temps local vers temps UTC
  Mm = str2int(strpart(date,15:16));
  S = str2int(strpart(date,18:19));
  //defines coordinates
  alpha = alpha2deg(alphaJ2000(1),alphaJ2000(2),alphaJ2000(3));
  delta = dec2deg(deltaJ2000(1),deltaJ2000(2),deltaJ2000(3));
  
  return airmassFromCoordinates(Y, M, D, H, Mm, S, alpha, delta);
}
func alpha2deg(h, m, s) {
        return dec2deg(h, m, s)*15;
}

func dec2deg(h, m, s) {
        return h+m/60.+s/3600.;
}
func airmassFromCoordinates(Y, M, D, H, Mm, S, alpha, delta, latitude, longitude)
/* DOCUMENT airmass(Y, M, D, H, Mm, S, alpha, delta, latitude, longitude)
     
   SEE ALSO:
 */
{
        LATITUDE = 28.7567;
        LONGITUDE = -17.8917;

        if (is_void(latitude)) latitude=LATITUDE;
        if (is_void(longitude)) longitude=LONGITUDE;

        alt = AAAltitude(Y, M, D, H, Mm, S, alpha, delta, latitude, longitude)(1);
        return (1/cos(pi/2.-alt*pi/180));
}
func AAAltitude(Y, M, D, H, Mm, S, alpha, delta, latitude, longitude) {
        LATITUDE = 28.7567;
        LONGITUDE = -17.8917;

        if (is_void(latitude)) latitude=LATITUDE;
        if (is_void(longitude)) longitude=LONGITUDE;


        lst = LST(Y, M, D, H, Mm, S, longitude) * 15.;
        H = lst - alpha;
        r = pi/180;
        sina = sin(latitude*r)*sin(delta*r)+cos(latitude*r)*cos(delta*r)*cos(H*r);
        alt = asin(sina) / r;
        return alt;
}
func LST(Y, M, D, H, Mm, S, longitude) {

        LONGITUDE = -17.8917;
        if (is_void(longitude)) longitude=LONGITUDE;

        JJ = JD(Y, M, D, H, Mm, S);
        D = JJ - 2451545.;
        GMST = 18.697374558 + 24.06570982441908 * D;
        lst = GMST+longitude/15.;
        lst%=24;
        return lst;
}

func JD(Y, M, D, H, Mm, S) {
        h=(H-12)/24. + Mm/24./60 + S/24./3600;
        J=(1461 * (Y + 4800 + (M - 14)/12))/4 +(367 * (M - 2 - 12 * ((M - 14)/12)))/12 - (3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075;
        return double(J) + h;
}


/*

 _   _  ____ ____   _    
| \ | |/ ___|  _ \ / \   
|  \| | |   | |_) / _ \  
| |\  | |___|  __/ ___ \ 
|_| \_|\____|_| /_/   \_\
                         
*/

func giveNCPADataDir(Dir){

  if(strpart(Dir,1:7)=="2013_09"){
    DirNCPA = "2013_09_15_onsky/";
    suffirbench = "19h57m41s";
  }
  else if(strpart(Dir,1:7)=="2013_07"){
    DirNCPA = "2013_07_22_onsky/";
    suffirbench = "17h41m59s";
  }

  else if( strpart(Dir,1:7)=="2012_10" ){
    DirNCPA = "2012_10_04_onsky/";
    suffirbench = "20h43m28s";
  }

  
  return [DirNCPA,suffirbench];
}

func getOTFncpa(nPixels,procDir,&SR_bench,&PSF_ncpa,disp=)
{

  // .... finding the directories where ncpa calibration are storaged
  tmp = giveNCPADataDir(procDir);
  DirNCPA = tmp(1);
  timencpa = tmp(2);
  dataDirRoot = dataDir +  DirNCPA;
  // .... loading, cleaning (background extraction, dead pixel processing) and cropping the best PSF got on bench
  PSF_ncpa = processCanaryPSF(timencpa,SR_bench,box=nPixels,disp=disp);
  OTF_ncpa = fft(roll(PSF_ncpa)).re;

  return OTF_ncpa;

}


func ncpaPhaseDiversity(procDir,nm,nlow,&psfNcpa,&a,verb=,disp=)
/* DOCUMENT psf2 = ncpaPhaseDiversity("2013_09_17_onsky/",20,36,psf,a)

 */
{
  N = cam.nPixelsCropped;
  OTF_ncpa = getOTFncpa(N,procDir,SR_bench,psfNcpa,disp=disp);


  require,"/home/omartin/CANARY/Algorithms/OPRA/opra.i";

  modes     = "zernike"; // 4.92e-3
  lambda    = cam.lambda;   // [m]
  pixsize   = cam.pixSize;     // [arcsec]
  nmodesmax = 100;

  opp = opra(psfNcpa(,,-:1:1),[0],lambda,pixsize,tel.diam,nmodes=nmodesmax,use_mode=modes, \
             noise=0.0,cobs=tel.obs,progressive=0,first_nofit_astig=1,fix_kern=0,fix_pix=1, \
  niter=10,fix_defoc=1,dpi=140,gui=0,nm=1);

/*
  //grabbing the PSF
  N = cam.nPixelsCropped;
  OTF_ncpa = getOTFncpa(N,procDir,SR_bench,psfNcpa,disp=disp);

  // .... GEOMETRY 

  //Defining the pupil
  x       = (indgen(N) -N/2)(,-:1:N);
  pixSize = tel.pixSize*(tel.nPixels/N);
  x       *= pixSize/(tel.diam/2.);
  y        = transpose(x);
  P        = circularPupFunction(x,y,tel.obs);

  // defining the polar corrdinates
  tmp   = polarFromCartesian(x,y);
  rho   = tmp(,,1);
  theta = tmp(,,2);

  
  //phase diversity

  input    = array(pointer,2);
  input(1) = &[nlow,N,pixSize];
  input(2) = &[P,rho,theta];
  a    = array(.1,nm);

  //Levenberg-Marquardt Iterative fitting
  res  = lmfit(psfModel,input,a,psfNcpa,eps=1e-3,verb=verb);

  return psfModel(input,a);
*/
}

func psfModel(input,a)
{

  nm      = numberof(a);

  if(is_pointer(input)){
    nlow    = int((*input(1))(1));
    N       = int((*input(1))(2));
    pixSize = (*input(1))(3);
    P       = (*input(2))(,,1);
    rho     = (*input(2))(,,2);
    theta   = (*input(2))(,,3);
  }else{
    N = cam.nPixelsCropped;
    nlow = 1;
    //Defining the pupil
    x       = (indgen(N) -N/2)(,-:1:N);
    pixSize = tel.pixSize*(tel.nPixels/N);
    x       *= pixSize/(tel.diam/2.);
    y        = transpose(x);
    P        = circularPupFunction(x,y,tel.obs);

    // defining the polar corrdinates
    tmp   = polarFromCartesian(x,y);
    rho   = tmp(,,1);
    theta = tmp(,,2);
  }
  
  //deriving the Zernike modes
  phi = array(0.,N,N);
  for(i=1;i<=nm;i++){
    phi += a(i) * computeZernikePolynomials(rho,theta,i+nlow) ;
  }

  //computing the psf
  psf  = abs(fft(roll(P*exp(1i*phi))))^2;

  //normalization
  k = 0.990707*(tel.diam^4 * (1. - tel.obs^2)^2 *pi*pi/16.)/pixSize^4;
  psf /= k;

  return roll(psf);
}

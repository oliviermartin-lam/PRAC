include, "mathUtils.i";
include, "contextUtils.i";
include, "constantUtils.i";
/*
 ____ _____ ____  _   _  ____ _____ _   _ ____  _____ ____  
/ ___|_   _|  _ \| | | |/ ___|_   _| | | |  _ \| ____/ ___| 
\___ \ | | | |_) | | | | |     | | | | | | |_) |  _| \___ \ 
 ___) || | |  _ <| |_| | |___  | | | |_| |  _ <| |___ ___) |
|____/ |_| |_| \_\\___/ \____| |_|  \___/|_| \_\_____|____/ 
                                                            
*/

struct srcFao_struct{
  double lambda;
  double magnitude;
  double zeroPoint;
  double zenith;
  double azimuth;
  double location(2);
};
struct telFao_struct{
  double D;
  double obs;
  double zenith;
  double aera;
  double aeraInPix;
};
struct atmFao_struct{
  double  lambda;
  int     nLayers;
  double  r0;
  double  altitude(NL);
  double  L0(NL);
  double  weight(NL);
  double  windSpeed(NL);
  double  windDir(NL);
};
struct rtcFao_struct{
  double loopGain;
  double Fe;
  double latency;
  string aoMode;
};

struct fao_struct{
  srcFao_struct src;
  telFao_struct tel;
  rtcFao_struct rtc;
  atmFao_struct atm;
  double  pitch;
  double  kc;
  int     resolution;
  int     nTimes;
  double  pixSize;
  double  wfsLambda;
  double  varNoise;
  pointer kt;
  pointer ktx;
  pointer kty;
  pointer k;
  pointer kx
  pointer ky;
  pointer idxPara;
  pointer idxOrth;
  pointer Rx;
  pointer Ry;
  pointer Gx;
  pointer Gy;
  pointer hrej;
  pointer hn;
  pointer ha;
  pointer spatialHrej;
  pointer spatialHn;
  pointer spatialHa;
  int     nTh;
  int     disp;
  pointer PSF;
  double  SR;
};

/*
 ___ _   _ ____ _____  _    _   _ _____ ___    _  _____ ___ ___  _   _ 
|_ _| \ | / ___|_   _|/ \  | \ | |_   _|_ _|  / \|_   _|_ _/ _ \| \ | |
 | ||  \| \___ \ | | / _ \ |  \| | | |  | |  / _ \ | |  | | | | |  \| |
 | || |\  |___) || |/ ___ \| |\  | | |  | | / ___ \| |  | | |_| | |\  |
|___|_| \_|____/ |_/_/   \_\_| \_| |_| |___/_/   \_\_| |___\___/|_| \_|
                                                                       
*/



func initFaoStruct(void)
/* DOCUMENT fao = initFaoStruct()

   Intialiazes the fao structures
   parameters.
   
 */
{
  NL                = 1;
  fao               = fao_struct();

  // SOURCES
  fao.src.lambda    = 1.6e-6;
  fao.src.magnitude = 12;
  fao.src.zeroPoint = 1e5;
  fao.src.zenith    = 3*arcsec2radian; //set the figure in arcsec
  fao.src.azimuth   = 0*arcsec2radian; 
  fao.src.location  = tan(fao.src.zenith) *[cos(fao.src.azimuth),sin(fao.src.azimuth)];
  
  // TELESCOPE
  fao.tel.D         = 4.2;
  fao.tel.obs       = .285;
  fao.tel.aera      = pi*fao.tel.D^2 * (1-fao.tel.obs^2);

  // ATMOSPHERE
  fao.atm.nLayers   = NL;
  fao.atm.r0        = .15;
  fao.atm.weight    = [1.];
  fao.atm.altitude  = [0.];
  fao.atm.L0        = [30.];
  fao.atm.windSpeed = [10.];
  fao.atm.windDir   = [0.];

  // RTC
  fao.rtc.Fe        = 150.;
  fao.rtc.loopGain  = .5;
  fao.rtc.latency   = 0/fao.rtc.Fe;
  fao.nTh           = 10;
  fao.rtc.aoMode    = "SCAO";

  // WFS NAD SYSTEM
  fao.pitch         = 0.6;
  fao.kc            = 1./2/fao.pitch;
  fao.resolution    = 81;
  fao.nTimes        = 3;
  fao.pixSize       = 2*fao.kc/fao.resolution;
  fao.wfsLambda     = .5e-6;
  fao.varNoise      = 0;

  // MISCALLANEOUS
  fao.disp          = 0;
  
  // RECONSTRUCTORS AND TRANSFER FUNCTIONS
  fao               = defineSpatialFrequency(fao);
  fao               = reconstructionFilter(fao);
  fao               = controller(fao);
  fao               = spatialTransferFunction(fao);

  
  //palette,"earth.gp";
 return fao;
}


func defineSpatialFrequency(fao)
/* DOCUMENT fao = defineSpatialFrequency(fao);

   Define the spatial Fourier frequency vectors and the index
   of un/correctable frequencies from fao.pitch.
   
 */
{
  N      = fao.resolution * fao.nTimes;
  // 1D frequency vector
  kx     = span(-1,1,N) * fao.kc * fao.nTimes;
  // Duplicating N times the first dimension to get a 2D map
  fao.kx = &kx(,-:1:N);
  fao.ky = &transpose(*fao.kx);
  fao.k  = &abs(*fao.kx,*fao.ky);
  // Frequency indexes
  fao.idxPara = &(abs(*fao.kx) <  fao.kc & abs(*fao.ky) <  fao.kc);
  fao.idxOrth = &(abs(*fao.kx) >= fao.kc | abs(*fao.ky) >= fao.kc);

  //Computing the truncated frequencies vector from -kc to kc only;
  N      = fao.resolution;
  ktx     = span(-1,1,N) * fao.kc;
  fao.ktx = &ktx(,-:1:N);
  fao.kty = &transpose(*fao.ktx);
  fao.kt  = &abs(*fao.ktx,*fao.kty);
    
  return fao;
}

/*
 ____  _____ ____ ___  _   _ ____ _____ ____  _   _  ____ _____ ___  ____  
|  _ \| ____/ ___/ _ \| \ | / ___|_   _|  _ \| | | |/ ___|_   _/ _ \|  _ \ 
| |_) |  _|| |  | | | |  \| \___ \ | | | |_) | | | | |     | || | | | |_) |
|  _ <| |__| |__| |_| | |\  |___) || | |  _ <| |_| | |___  | || |_| |  _ < 
|_| \_\_____\____\___/|_| \_|____/ |_| |_| \_\\___/ \____| |_| \___/|_| \_\
                                                                           
*/

func reconstructionFilter(fao)
/* DOCUMENT fao= reconstructionFilter(fao)

   Computes the spatial filtering matrix operated by the Shack-Hartmann
   WFS.
 */
{

  //Grabbing parameters
  d   = fao.pitch
  Ts  = 1./fao.rtc.Fe;
  td  = fao.rtc.latency;
  kx  = *fao.ktx;
  ky  = *fao.kty;

  // Shack Hartmann WFS operator (Gradient)
  Gx = 1i*2*pi*kx*d*sinc(d*kx)*sinc(d*ky)*exp(1i*pi*d*(kx+ky));
  Gy = 1i*2*pi*ky*d*sinc(d*kx)*sinc(d*ky)*exp(1i*pi*d*(kx+ky));
  G  = Gx^2 + Gy^2;

  // Computing the phase reconstructor : Wiener filtering with gamma = 1
  aPSD       = atmospherePSD(fao); // atmosphere PSD
  ki         = fao.resolution*(fao.nTimes-1)/2 + 1;
  kf         = fao.resolution*(fao.nTimes+1)/2;
  aPSD       = aPSD(ki:kf,ki:kf);
  nPSD       = fao.varNoise/(2*fao.kc)^2; // noise PSD

  // Reconstructor derivation
  Rx         = conj(Gx)*aPSD/(G*aPSD + nPSD + 1e-15);
  Ry         = conj(Gy)*aPSD/(G*aPSD + nPSD + 1e-15);
  N          = dimsof(Rx)(0);
  Rx(N/2+1,N/2+1) = 0;
  Ry(N/2+1,N/2+1) = 0;

  //Storaging into fao struct
  fao.Gx = &Gx;
  fao.Gy = &Gy;
  fao.Rx = &Rx;
  fao.Ry = &Ry;
  
  return fao;
}

/*
  ____ ___  _   _ _____ ____   ___  _     _     _____ ____  
 / ___/ _ \| \ | |_   _|  _ \ / _ \| |   | |   | ____|  _ \ 
| |  | | | |  \| | | | | |_) | | | | |   | |   |  _| | |_) |
| |__| |_| | |\  | | | |  _ <| |_| | |___| |___| |___|  _ < 
 \____\___/|_| \_| |_| |_| \_\\___/|_____|_____|_____|_| \_\
                                                            

*/

func controller(fao)
/* DOCUMENT fao = controller(fao)

 */
{
  // Parameters
  Te    = 1./fao.rtc.Fe;
  delta = fao.rtc.latency/Te;
  delay = 1 + delta;
  g     = fao.rtc.loopGain;
  //define the Laplace variable s and the discrete variable z
  //nu = spanl(-3,log10(.5/Te),500);
  nu = (span(0,.5/Te,500))(2:);
  s  = 2i*pi*nu;
  z  = exp(s*Te);

  // WFS+DAC  transfer function
  hwfs =  z^(-delay);
  //defines the transfer functions
  if(fao.rtc.aoMode == "MOAO"){
    hInt = g*z/(z-1+g);
    hrej = 1. - hInt*hwfs;
    hn   = hInt*hwfs/z;
    ha   = hInt*hwfs;
  }else{
    hInt = g*z/(z - 1.);
    hrej = 1./(1. + hInt*hwfs);
    hn   = hInt*hwfs/z/(1. + hInt*hwfs);
    ha   = hInt*hwfs/(1. + hInt*hwfs);
  }
 
  fao.hrej = &hrej;
  fao.hn   = &hn;
  fao.ha   = &ha;
  
  //Display
  if(fao.disp){
    bodi = bode(hInt);
    bodr =  bode(hrej);
    bodn =  bode(hn);
    boda =  bode(ha);

    winkill,0;
    window,0,dpi=90,style="aanda.gs";fma;limits;logxy,1,1;gridxy,1,1;
    plg,bodr(,1),nu,marks=0;
    plg,bodn(,1),nu,marks=0,color="red";
    plg,boda(,1),nu,marks=0,color="blue";
    plg,bodi(,1),nu,marks=0,color="green";
    xytitles,"Temporal frequency [Hz]","Transfer functions [dB]";
    plt,"Signal rejection",0.2,10.,tosys=1;
    plt,"Noise rejection",0.2,5.,tosys=1,color="red";
    plt,"Aliasing rejection",0.2,2,tosys=1,color="blue";
    plt,"Integrator",0.2,1,tosys=1,color="green";

    winkill,1;
    window,1,dpi=90,style="aanda.gs";fma;limits;logxy,1,0;gridxy,1,1;
    plg,bodr(,2),nu,marks=0;
    plg,bodn(,2),nu,marks=0,color="red";
    plg,boda(,2),nu,marks=0,color="blue";
    xytitles,"Temporal frequency [Hz]","Phase [Deg]";
    plt,"Signal rejection",0.2,70.,tosys=1;
    plt,"Noise rejection",0.2,50.,tosys=1,color="red";
    plt,"Aliasing rejection",0.2,30,tosys=1,color="blue";
    
  }
  

  return fao;
}

func spatialTransferFunction(fao)
/* DOCUMENT fao = spatialTransferFunction(fao)

   Computes the equivalent spatial transfer function from the temporal
   transfer function averaged over layers and directions.
 */
{
// Parameters
  Te    = 1./fao.rtc.Fe;
  delta = fao.rtc.latency/Te;
  delay = 1 + delta;
  g     = fao.rtc.loopGain;
  
  //Equivalent spatial transfer function
  kx = *fao.ktx;
  ky = *fao.kty;
  vx = fao.atm.windSpeed*cos(fao.atm.windDir*pi/180.);
  vy = fao.atm.windSpeed*sin(fao.atm.windDir*pi/180.);
  Hrej = Hn = Ha = 0*kx;
  
  // Loop on discrete turbulent layers
  for(l=1;l<=fao.atm.nLayers;l++){
    //Loop on direction angle
    for(t=1;t<=fao.nTh;t++){
      theta = 2*pi*(t-1)/fao.nTh;
      Nu    = abs(vx(l)*kx*cos(theta),vy(l)*ky*sin(theta));
      s     = 2i*pi*Nu;
      z     = exp(s*Te);
      hwfs =  z^(-delay);
      
      if(fao.rtc.aoMode == "MOAO"){
        hInt  = g*(delta + (1.-delta)*z)/(z*(z-1+g));
        Hrej += fao.atm.weight(l)*(1. - hInt*hwfs)/fao.nTh;
        Hn   += fao.atm.weight(l)*(hInt*hwfs/z)/fao.nTh;
        Ha   += fao.atm.weight(l)*(hInt*hwfs)/fao.nTh;
      }else{
        hInt  = g*(delta + (1.-delta)*z)/(z*(z-1));
        Hrej += fao.atm.weight(l)*(1./(1 + hInt*hwfs))/fao.nTh;
        Hn   += fao.atm.weight(l)*(hInt*hwfs/z/(1. + hInt*hwfs))/fao.nTh;
        Ha   += fao.atm.weight(l)*(hInt*hwfs/(1. + hInt*hwfs))/fao.nTh;
      }
    }
  }
  //Storaging
  fao.spatialHrej = &Hrej;
  fao.spatialHn   = &Hn;
  fao.spatialHa   = &Ha;
  return fao;
}


/*
 ___ _____ _____ ______  ____  _____ 
 / _ \_   _|  ___/ /  _ \/ ___||  ___|
| | | || | | |_ / /| |_) \___ \| |_   
| |_| || | |  _/ / |  __/ ___) |  _|  
 \___/ |_| |_|/_/  |_|   |____/|_|    
                                      
                  
*/
func aoPSF(fao,ol=)
/* DOCUMENT PSF = aoPSF(fao,ol=)

   Computes and displays the long-exposure time PSF delivered
   by the AO system.
*/
{
  // Grabbing the OTF
  OTF = roll(aoOTF(fao,ol=ol));
  //Computing the PSF
  fao.PSF   = &(roll(fft(OTF).re)/fao.pixSize);
  *fao.PSF /= fao.tel.aera;
  //Strehl ratio computation
  fao.SR  = 100*sum(OTF)/sum(telescopeOTF(fao));
  //Display
  window,2;fma;limits;palette,"earth.gp";
  a = fao.pixSize * fao.resolution * fao.nTimes;
  pli,*fao.PSF,-a,-a,a,a;
  xytitles,"Arcsec","Arcsec";
  pltitle,"Estimated PSF\n Strehl = " + var2str(arrondi(fao.SR,1)) + "%";

  //Shows the DM cut-off frequency limit
  kc = fao.kc * fao.src.lambda * radian2arcsec;
  plg,[-kc,-kc],[-kc,kc] ,marks=0,type=2,color="white";
  plg,[-kc,kc] ,[kc,kc]  ,marks=0,type=2,color="white";
  plg,[kc,kc]  ,[kc,-kc] ,marks=0,type=2,color="white";
  plg,[kc,-kc] ,[-kc,-kc],marks=0,type=2,color="white";

}
func aoOTF(fao,ol=)
/* DOCUMENT OTF = aoOTF(fao,ol=)

   Computes the OTF from the residual phase.
*/
{

  // Grabbing the PSD from the AO system
  PSD = aoPSD(fao,ol=ol);
  //Computes the OTF
  sysOTF = PSD2OTF(PSD,fao.pixSize);
  // Grabbing the telescope OTF
  telOTF = telescopeOTF(fao);

  return telOTF * sysOTF;
}

/*
 _____ _____ _     _____ ____   ____ ___  ____  _____ 
|_   _| ____| |   | ____/ ___| / ___/ _ \|  _ \| ____|
  | | |  _| | |   |  _| \___ \| |  | | | | |_) |  _|  
  | | | |___| |___| |___ ___) | |__| |_| |  __/| |___ 
  |_| |_____|_____|_____|____/ \____\___/|_|   |_____|
                                                      
*/


func telescopeOTF(fao)
/* DOCUMENT OTF = telescopeOTF(fao)

   Returns the perfect Telescope OTF with at a resolution
   equal to fao.resolution * fao.nTimes.
*/
{
  // Defining geometry
  N     = fao.resolution * fao.nTimes;
  x     = (indgen(N) -N/2)(,-:1:N);
  x    *= fao.pixSize/(fao.tel.D/2.);
  y     = transpose(x);
  // Pupil definition
  r = abs(x,y);
  P = r<=1.0 & r>fao.tel.obs;
  // FFT autocorrelation
  OTF = fft(abs(fft(P))^2).re / dimsof(P)(0)^2;

  // Normalize the OTF to get a psf with PSF(0)=1.00 when diffraction-limited
  OTF /= sum(OTF);

  return roll(OTF);

}
/*
 ____  ____  ____  
|  _ \/ ___||  _ \ 
| |_) \___ \| | | |
|  __/ ___) | |_| |
|_|   |____/|____/ 
                   
*/
func aoPSD(fao,ol=)
/* DOCUMENT PSD = aoPSD(fao,ol=)

   Returns the PSD from the AO system.
*/
{
  //Initialization
  N = fao.resolution * fao.nTimes;
  PSD = array(0.,N,N);
  
  //Adding low frequencies PSD
  ki = (N - fao.resolution)/2 + 1;
  kf = (N + fao.resolution)/2;
  PSD(ki:kf,ki:kf) += 0*servoLagPSD(fao) + 0*aliasingPSD(fao,ol=ol) + 0*noisePSD(fao,ol=ol);
  
  //Adding fitting PSD
  PSD += fittingPSD(fao);

  return PSD;
}

func atmospherePSD(fao)
/* DOCUMENT PSD = atmospherePSD(fao);

   Computes the atmospheric phase PSD using the Wiener
   spectrum expressed at any layers and summed up
   over altitude.

   It invokes the atmosphere properties of the fao struct:
   fao.r0,fao.weight, fao.altitude and fao.L0 computed at
   fao.lambda.
 */
{
  N      = fao.resolution * fao.nTimes;
  PSD    = array(0.,N,N);
  nl     = fao.atm.nLayers;
  k      = *fao.k;
  r0     = fao.atm.r0 * (fao.src.lambda/fao.wfsLambda)^(6./5);
  L0     = fao.atm.L0;
  weight = fao.atm.weight;

  cte = (24*gamma_R(6./5)/5)^(5./6)*(gamma_R(11./6)^2/(2.*pi^(11./3)));

  // Summation over altitude
  for(l=1;l<=nl;l++){
    PSD += cte * r0^(-5/3.) * weight(l) * (k^2. + 1./L0(l)^2.)^(-11./6);
  }
         
           
  return PSD;
}

func fittingPSD(fao)
/* DOCUMENT PSD = fittingPSD(fao)

   Computes the fitting PSD from the atmosphere PSD
   by removing frequencies lower than fao.kc. The DM footprint
   in the Fourier domain is assumed to be a square.
*/
{
  PSD = atmospherePSD(fao);
  PSD = PSD * (*fao.idxOrth);
  //PSD = PSD * pistonRemoval(fao.tel.D,fao.pitch,*fao.k,0,0);
  return PSD;
}

func servoLagPSD(fao)
/* DOCUMENT PSD = servoLagPSD(fao);

   Returns the servo-Lag PSD from the latency error introduced
   by the AO system.
   
 */
{

  //Spatial transfer function
  H = *fao.spatialHrej;
  
  //Atmospheric phase PSD
  w         = where(*fao.idxPara);
  atmPSD    = 0*H;
  atmPSD(*) = (atmospherePSD(fao))(w);
     
  //AO transfer function
  PSD = H*atmPSD;

  // Removing Piston and high spatial frequencies
  PSD = PSD * pistonRemoval(fao.tel.D,fao.pitch,*fao.kt,0,0);

  return abs(PSD);
}


func noisePSD(fao,ol=)
/* DOCUMENT PSD = noisePSD(fao,varNoise)

   Returns the noise spatial PSD.
*/
{
  if(is_void(ol))
    ol = 0;
  
  kx   = *fao.ktx;
  ky   = *fao.kty;
  d    = fao.pitch;
  varn = fao.varNoise/fao.kc;

  //Noise PSD
  PSD  = varn/(2*fao.kc)^2 *(abs(*fao.Rx)^2 + abs(*fao.Ry)^2);
  
  if(!ol){
    //Filtering
    PSD *=  abs(*fao.spatialHn)^2;
  }
    
  PSD = PSD * pistonRemoval(fao.tel.D,fao.pitch,*fao.kt,0,0);

  return PSD;
}

func aliasingPSD(fao,ol=)
/* DOCUMENT PSD = aliasingPSD(fao)

   Returns the aliasing spatial PSD;
*/
{
  //Grabbing parameters
  fr0    = fao.atm.weight;
  T      = 1./fao.rtc.Fe;
  td     = fao.rtc.latency;
  D      = fao.tel.D;
  d      = fao.pitch;
  w      = 2i*pi*d;
  klim   = fao.nTimes;
  kx     = *fao.ktx;
  ky     = *fao.kty;
  k      = *fao.kt;

  if(ol)
    ha = 1;
  else
    ha = *fao.spatialHa;
  
  PSD    = 0*kx;
  
  // Effective outer scale calculation
  L0     = sum(fao.atm.weight * fao.atm.L0^(-5/3.))^(-3/5.);
  f0     = 1./L0;
  //Wind speed distribution
  vx     = fao.atm.windSpeed*cos(fao.atm.windDir*pi/180.);
  vy     = fao.atm.windSpeed*sin(fao.atm.windDir*pi/180.);

  // Double on x/y aliased spatial frequencies
  for(m=-klim;m<=klim;m++){
    for(n=-klim;n<=klim;n++){

      if(m!=0 || n!=0){
        //Atmospheric phase PSD
        km   = kx - m/d;
        kn   = ky - n/d;
        W_mn = (km^2 + kn^2 + f0^2)^(-11/6.);
        W_mn = W_mn *  pistonRemoval(D,d,k,m,n);
        //Spatial reconstructor
        Q = ((*fao.Rx)*w*km + (*fao.Ry)*w*kn) * (sinc(d*km)*sinc(d*kn));

        //Loop on layers
        avr = 0;
        for(l=1;l<=fao.atm.nLayers;l++){
          delta = exp(2i*pi*km*vx(l)*td)*exp(2i*pi*kn*vy(l)*td);
          avr  += fr0(l)* sinc(km*vx(l)*T) * sinc(kn*vy(l)*T) * delta;
        }
        q    = Q* avr*ha;
        PSD += W_mn * (q*conj(q));
      }
    }
  }
  
  // Removing Piston
  r0  = fao.atm.r0 * (fao.src.lambda/fao.wfsLambda)^(6/5.);
  PSD*= 0.0229 * r0^(-5/3.);
  PSD = PSD * pistonRemoval(fao.tel.D,fao.pitch,k,0,0);
 
  return PSD.re;
}

func anisoplanatismPSD(fao)
/* DOCUMENT PSD =  anisoplanatismPSD(fao)

   Returns the atmosphere PSD multiplied by the
   anisoplanatism filter;
 */
{
  zLayer = fao.atm.altitude;
  fr0    = fao.atm.weight;
  kx     = *fao.kx;
  ky     = *fao.ky;
  k      = *fao.k;
  alpha  = fao.src.location;
  N      = dimsof(kx)(0);
  A      = array(0.,N,N);

  for(l=1;l<=fao.atm.nLayers;l++){
    red  = 2*pi*zLayer(l)* ( kx*alpha(1) + ky*alpha(2) );
    A   += fr0(l)*( 1 - cos(red) );
  }

  PSD = A * atmospherePSD(fao) * pistonRemoval(fao.tel.D,fao.pitch,k,0,0);
  
  return PSD;
}

func theoriticalErrorBreakdown(fao)
/* DOCUMENT out = theoriticalErrorBreakdown(fao);

   Computes analytically WF variance of each individual system error.
*/
{
  rad2nm = 1e9*fao.src.lambda/2/pi;

  // Fitting
  r0 = fao.r0*(fao.src.lambda/fao.wfsLambda)^(6/5.);
  wfeFit = sqrt(0.23*(fao.pitch/r0)^(5/3.))*rad2nm;

  //Aliasing
  wfeAlias = sqrt(0.074*(fao.pitch/r0)^(5/3.))*rad2nm;

  //Noise
  wfeNoise = sqrt(fao.varNoise)*rad2nm;

  //Servo-lag
  //wfeServo = 

  return [wfeFit,wfeAlias,wfeNoise];
}

//////////////////////////////// ROUTINES ///////////////////////////////////////////

func PSD2COV(PSD,pixSize)
/* DOCUMENT COV = PSD2COV(PSD)

   Computes the spatial covariance map from the spatial PSD.
*/
{
  return roll(fft(roll(PSD)))*pixSize^2;
}
func COV2SF(COV)
/* DOCUMENT SF = COV2SF(COV)

   Computes the spatial phase structure function  map from thecovariance
   map COV
*/
{
  N = dimsof(COV)(0);
  return 2*(COV(N/2+1,N/2+1) - COV);
}
func SF2OTF(SF)
/* DOCUMENT OTF = SF2OTF(SF)

   Computes the OTF from the spatial phase structure function.
*/
{
  return exp(-.5*SF);
}
func PSD2OTF(PSD,pixSize)
/* DOCUMENT OTF = PSD2OTF(PSD,pixSize)

   Computes the OTF from the sptial PSD.
*/
{
  COV = PSD2COV(PSD,pixSize);
  SF  = COV2SF(COV);
  return SF2OTF(SF).re;
  
}

func pistonRemoval(D,d,k,m,n)
/* DOCUMENT PR = pistonRemoval(D,d,f,mi,ni)

   Returns the piston filter on spatial PSD
*/
{
  
  Bx = pi*D*(k-m/d); 
  By = pi*D*(k-n/d);
  B  = abs(Bx,By);
  F  = 2*bessj(1,B)/B;

 return  1 - abs(F)^2;
}

func wfe2sr(wfe,lambda)
/* DOCUMENT SR = wfe2sr(wfe,lambda);

   Determines the Strehl ratio at lambda from wfe in nm
   using the Marechal approximation
*/
{
  return exp(-(2*pi*wfe*1e-9/lambda)^2)
}

func sr2wfe(SR,lambda)
/* DOCUMENT wfe = sr2wfe(SR,lambda);

   Determines the Wave Front error in nm from the
   Strehl ratio at lambda  using the Marechal approximation
*/
{
  return sqrt(-log(SR))*lambda*1e9/2/pi;
}

func clr(void)
/* DOCUMENT clr;

   Do "fma" followed by "limits".
*/
{
  fma;limits;
}

func bode(h)
/* DOCUMENT [MAG,ARG] = bode(h)

   Returns the attenuation in dB and the
   phase in degree of the transfer function h.
*/
{

  // Magnitude
  MAG = 20*log10(abs(h));

  //PHASE
  x   = h.re;
  y   = h.im;
  ARG = 0*MAG;

  w0 = where(x!=0);
  if(is_array(w0))
    ARG(w0) = atan(y(w0)/x(w0));
    

  return[MAG,ARG*180/pi];
}

func filteringNoiseFactor(g,delay,obsmode)
/* DOCUMENT a = filteringNoiseFactor(g,delay,obsmode)

   Returns the propagation coefficient on a white noise
   propagating through a 1rst order low-pass filter
*/
  
{
  if(g==0) return 0;
  
  if(obsmode =="MOAO"){
    return g*(1-2*g*delay*(1-delay))/(2.-g);
  }else if(obsmode == "SCAO"){
    return g*g/(g*(1-delay)*(2 - g*(1-delay) + 2*g*delay*(1-g*(1-delay))^2/(1+g*delay) - (g*delay)^2 ));
  }
}

//fao = initFaoStruct();

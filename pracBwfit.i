func computeOTFbandwidth(geometry,obsMode,verb=)
/* DOCUMENT
 */
{
  // .... PSD computing using the Wiener expression
  r0tot = atm.r0*(cam.lambda/atm.lambda)^1.2;
  k = computeSpatialFreqRad(tel.nPixels, tel.fourierPixSize, kx, ky);
  Wiener = computeWienerSpectrum(k, tel.nPixels, r0tot, atm.L0);
  
  // .... Defining the frequency aera correctable by the DM
  msk = defineDmFrequencyArea(k, kx, ky, geometry, dm.pitch );
  ndm    = where( msk );
  northo = where( !msk );

  // .... Multipling by the system Z-transfer function at each altitude
  PSD_bw = computePSDbandwidth(k,rtc.Fe,rtc.delay,rtc.loopGain,rtc.BP,Wiener,northo,mode=obsMode,verb=verb);

  // .... Retrieving the phase structure function
  N = tel.nPixels;
  PSD_bw(N/2+1,N/2+1) = 0.0;
  PSD_bw(N/2+1,N/2+1) = -sum(PSD_bw);  
  DPHI_bw = 2*abs(fft(PSD_bw)) *  (tel.fourierPixSize^2);
  
  OTF_bw = exp(-0.5*DPHI_bw);

  return OTF_bw;
}
func computePSDbandwidth(k,Fe,tret,gain,BP,Wiener,northo,mode=,verb=)
/* DOCUMENT Wbp = computeBpSpectrum(k,V,dir,Fe,tret,gain,Wiener,northo)
     
   Bandwidth error. We use the transfer function of the system,
   defined by a sampling frequency and delay, assuming a pure
   close-loop integrator. The attenuation of the transfer function
   Hcor(nu) is applied on all spatial frequency k=nu/V (units equation
   is m^-1 = s^-1 / (m/s) ).

   SEE ALSO:
 */
{


  PSD_bw   = 0*Wiener;
  r0tot = atm.r0*(cam.lambda/atm.lambda)^1.2;
  cnh   = atm.cnh;
  l0    = atm.l0h;
  vh    = 3*atm.vh;
  //Compensation of the model unaccuracy of the spatial covariance
  cnh *= r0tot^(-5/3.)/sum(cnh);

  for(l = 1;l<=atm.nLayers;l++){
    nu_l = k * vh(l)/sqrt(2);
    if(mode == "MOAO"){
      hcor = hcorMoao(nu_l, Fe, tret, gain, BP);
    }else if(mode == "SCAO"){
      hcor = hcorScao(nu_l, Fe, tret, gain, BP);
    }
    PSD_bw +=  abs(hcor) * 0.023 * cnh(l) * (k^2 + 1./l0(l)^2.)^(-11/6.);
  }

  N = tel.nPixels;
  PSD_bw(N/2+1,N/2+1) = 0;
  PSD_bw(northo) = 0.00;

  return PSD_bw;
}


func computeOTFfitting(geometry,verb=)
/* DOCUMENT
 */
{
  r0tot = atm.r0*(cam.lambda/atm.lambda)^1.2;
  N = tel.nPixels;

  if( geometry == "circle" ) {
    //computing r0 at cam.lambda
    x = roll( ( indgen(N)-(N/2+1) ) * tel.pixSize );
    x = x(,-:1:N);                
    y = transpose(x);             
    r = abs(x,y);         
    DPHI_fit = dphi_highpass(r,tel.pitch) * r0tot^(-5./3);

  }else if(geometry == "square") {
    // initializing with the full frequency PSD
    k = computeSpatialFreqRad(tel.nPixels, tel.fourierPixSize, kx, ky);
    msk = defineDmFrequencyArea(k, kx, ky, geometry, dm.pitch );
    //defining the 0/1 square mask
    ndm    = where( msk );
    PSD_fit = computeWienerSpectrum(k, tel.nPixels, r0tot, atm.L0);
    PSD_fit(ndm) = 0.0;
    PSD_fit(N/2+1,N/2+1) = 0.0;
    PSD_fit(N/2+1,N/2+1) = -sum(PSD_fit);  
    DPHI_fit = 2*abs(fft(PSD_fit)) *  (tel.fourierPixSize^2);
  }else{

    //foutier frequencies
    k = computeSpatialFreqRad(tel.nPixels, tel.fourierPixSize, kx, ky);
    
    //compute the central influence function
    x = roll( ( indgen(N)-(N/2+1) ) * tel.pixSize );
    x = x(,-:1:N);                
    y = transpose(x);
    iF = funcInflu(x,y,dm.x0);

    //defining the influence transfer function
    tmp = readfits("fitsFiles/influenceFilter.fits",err=1);

    if(is_void(tmp)){
      H0 = fft(iF);
      Hif = array(0.,N,N);
      for(i=-dm.nActu;i<=dm.nActu;i++){
        for(j=-dm.nActu;j<=dm.nActu;j++){
          Hif += roll(H0,[int(2*i/tel.pitch),int(2*j/tel.pitch)]);
        }
      }
      writefits,"fitsFiles/influenceFilter.fits",[Hif.re,Hif.im];
    }else{
      Hif = tmp(,,1) + 1i*tmp(,,2);
    }
    
    //defining the fourier filter
    Hif    /= max(abs(Hif));
    Hfilter = abs(1 - Hif)^2;
    //computing the fitting PSD
    PSD_fit = computeWienerSpectrum(k, tel.nPixels, r0tot, atm.L0);
    PSD_fit *= roll(Hfilter);
    PSD_fit(N/2+1,N/2+1) = 0.0;
    PSD_fit(N/2+1,N/2+1) = -sum(PSD_fit);  
    //computing the structure function
    DPHI_fit = 2*abs(fft(PSD_fit)) *  (tel.fourierPixSize^2);
     
  }
  

  OTF_fit = exp(-0.5*DPHI_fit);
 
  return OTF_fit;
}

func computeSpatialFreqRad(N, uk, &kx, &ky)
{
  kx = ( indgen(N)-(N/2+1) ) * uk;   // on genere la composante X du vecteur freq spatiale k
  kx = kx(,-:1:N);                   // en 2 dims ...
  ky = transpose(kx);                // et pareil pour Y
  k = abs(kx,ky);
  k(N/2+1,N/2+1) = uk;           // valeur non-nulle pour ne pas calculer 0.00^(-11/3)=infini
  return k;
}

func defineDmFrequencyArea(k, kx, ky, sgeom, pitch )
{
  fc = 0.5/pitch;
  
  if( sgeom=="circ" )
    msk = k < fc;
  else
    msk = (kx < fc) & (kx > -fc) & (ky < fc) & (ky > -fc);
  return msk;
}

func computeWienerSpectrum(k, N, r0tot, outerScale)
{
  Wiener = 0.023 * r0tot^(-5./3) * (k*k + 1./outerScale^2.)^(-11./6);
  // frequence 0 mise a 0, vu que de ttes facons elle est `fausse'
  Wiener(N/2+1,N/2+1)=0;         
  return Wiener;
}


/*
 _____ ____      _    _   _ ____  _____ _____ ____  
|_   _|  _ \    / \  | \ | / ___||  ___| ____|  _ \ 
  | | | |_) |  / _ \ |  \| \___ \| |_  |  _| | |_) |
  | | |  _ <  / ___ \| |\  |___) |  _| | |___|  _ < 
  |_| |_| \_\/_/   \_\_| \_|____/|_|   |_____|_| \_\
                                                    
 _____ _   _ _   _  ____ _____ ___ ___  _   _ ____  
|  ___| | | | \ | |/ ___|_   _|_ _/ _ \| \ | / ___| 
| |_  | | | |  \| | |     | |  | | | | |  \| \___ \ 
|  _| | |_| | |\  | |___  | |  | | |_| | |\  |___) |
|_|    \___/|_| \_|\____| |_| |___\___/|_| \_|____/ 
                                                    
*/

func filteringNoiseFactor(g,delay,obsmode)
{
  if(g==0) return 0;
  if(obsmode =="MOAO"){
    return g*(1-2*g*delay*(1-delay))/(2.-g);
  }else if(obsmode == "SCAO"){
    return g*g/(g*(1-delay)*(2 - g*(1-delay) + 2*g*delay*(1-g*(1-delay))^2/(1+g*delay) - (g*delay)^2 ));
  }
}

func hsysScao(freq,Fe,tret,G)
/* DOCUMENT hsys = hsysScao(freq,Fe,tret,G)

   Returns the transfer function between the modes of the on-axis
   residual phase and the  on-axis phase. It depends on the frequencies domain freq, the
   sampling frequency freq (Hz), the gain of the loop and the delay
   tret (seconds) between the end of the integration of the WFS and
   the applying of the command.

   Olivier Martin.
    
   SEE ALSO: hwfs,hbfScao,hcorScao,hdm
 */
{
  
  Te = 1./Fe;
  //fractional delay in number of frames.
  delta = tret*Fe;
  //define the discrete z variable
  p = 2i*pi*freq + 1e-12;
  z = exp(p*Te);
  //defines the transfer function
  hsys = G*(delta + (1.-delta)*z)/(z*(z-1));

  return hsys;
}

func hbfScao(freq,Fe,tret,G,BP)
/* DOCUMENT hsys = hbfScao(freq,Fe,tret,G,BP)

   Returns the transfer function between the modes of the perpendicular
   phase and the on-axis phase. It depends on the frequencies domain freq, the
   sampling frequency freq (Hz), the gain of the loop and the delay
   tret (seconds) between the end of the integration of the WFS and
   the applying of the command.

   Olivier Martin.
    
   SEE ALSO: hwfs,hsysScao,hcorScao,hdm
 */
{
  hbo = hsysScao(freq,Fe,tret,G) * hdm(freq,BP) * hwfs(freq,Fe);
  return hbo/(1. + hbo);
}

func hcorScao(freq,Fe,tret,G,BP)
/* DOCUMENT   hcor(freq,Fe,tret,G,BP)

   
   The option an=1 sets the integrator to "analog". Doing this, an
   extra 1/2 frame delay is added compared to case of the numeric
   integrator an=0.

   <tret> is the delay expressed as a *time in seconds*, between the
   end of the integration and the start of the command.
   
*/
{
  Te=1./Fe;
  p = 1i*2*pi*freq + 1e-12;
    
  hbo = hsysScao(freq,Fe,tret,G) * hdm(freq,BP) * hwfs(freq,Fe);

  hcor = 1./(1. + hbo);
  
  return hcor;
}


func hdm(freq,BP)
/* DOCUMENT hmir =  hdm(freq,BP)

   Retuns the transfer function of the mirror. It depends on the
   frequencies domain freq, the sampling frequency freq (Hz) and the
   bandwidth of the mirror BP (Hz).

   Olivier Martin.
   
   SEE ALSO: hsysMoao,hwfs,hboMoao,hcorMoao
*/
{
  return 1./(1. + 1i*freq/BP);
}

func hwfs(freq,Fe)
/* DOCUMENT haso = hwfs(freq,Fe)

   Retuns the transfer function of the WFS. It depends on the
   frequencies domain freq, the sampling frequency Freq (Hz), the gain
   of the loop and the delay tret (seconds) between the end of the
   integration of the WFS and the applying of the command.

   Olivier Martin.

   SEE ALSO: hsysMoao,hboMoao,hcorMoao,hdm
 */
{

  Te = 1./Fe;
  //define the discrete z variable
  p = 1i*2*pi*freq + 1e-12;
  z = exp(-p*Te);
  //defines the transfer function
  haso =  (1. - z)/(p*Te);
  haso*=haso; //taking into account both CCD and DAC delay

  return haso;
}

func hsysMoao(freq,Fe,tret,G)
/* DOCUMENT hsys = hsysMoao(freq,Fe,tret,G)

   Returns the transfer function between the modes of the on-axis
   residual phase and the estimated on-axis phase from the tomographic
   reconstruction. It depends on the frequencies domain freq, the
   sampling frequency freq (Hz), the gain of the loop and the delay
   tret (seconds) between the end of the integration of the WFS and
   the applying of the command.

   Olivier Martin.
    
   SEE ALSO: hwfs,hboMoao,hcorMoao,hdm
 */
{
  
  Te = 1./Fe;
  delta = tret*Fe;//fractionnal delay in number of frames.
  //define the discrete z variable
  p = 1i*2*pi*freq + 1e-12;
  z = exp(-p*Te);
  //defines the transfer function
  hsys = G*(delta + (1.-delta)*z)/(z*(z-1+G));

  return hsys;
}


func hboMoao(freq,Fe,tret,G,BP)
/*
  Retuns the transfer function of the MOAO system in open-loop. It
  depends on the frequencies domain freq, the sampling frequency freq
  (Hz), the gain of the loop, the bandwith of the mirror BP (Hz) and
  the delay tret (seconds) between the end of the integration of the
  WFS and the applying of the command. We have:

  hboMoao = hsysMoao*haso*hdm*hdac with hdac = haso

  Olivier Martin.
   
  SEE ALSO: hsysMoao,hwfs,hcorMoao,hdm
 */
{
  
  return  hsysMoao(freq,Fe,tret,G)*hdm(freq,BP) *hwfs(freq,Fe);
}

func hcorMoao(freq,Fe,tret,G,BP)
/*
  Retuns the transfer function of the MOAO system in engaged-loop. It
  depends on the frequencies domain freq, the sampling frequency freq
  (Hz), the gain of the loop, the bandwidth of the mirror BP (Hz) and
  the delay tret (seconds) between the end of the integration of the
  WFS and the applying of the command. We have :

  beta: interest directions (TS for Canary)
  alpha: guiding directions
  
  e_para(beta) = hcorMoao*R*a_para(alpha) - hboMoao*R*a_ortho(alpha)
                - hsysMoao*R*n(alpha)+ [a_para(beta)-R*a_para(alpha)]
  
  e_para(beta): modes of the residual parallel phase in the beta direction
  a_para(alpha): modes of the atmospherical parallel phase in the alpha direction
  a_ortho(alpha): modes of the atmospherical orthogonal phase in the alpha direction
  n(alpha): modes of the noise measured by the off-axis WFS
  R: tomographic reconstructor

  Olivier Martin.
  
  SEE ALSO: hsysMoao,hwfs,hboMoao,hmirror
 */
{
  return 1. - hboMoao(freq,Fe,tret,G,BP);
}


func psd2rms( W )
{
  fact = tel.fourierPixSize * cam.lambda /2/pi * 1e9;
  return sqrt(sum(W)) * fact;
}

func dphi2rms(dphi)
{
  fact =  cam.lambda /2/pi * 1e9;
  return sqrt(sum(abs(fft(-0.5*cam.dPixSizeOnLambda*dphi)/(numberof(dphi)-1.))))*fact;
}






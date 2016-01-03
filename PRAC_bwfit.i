/*
 ____  _   _    _    ____  _____ 
|  _ \| | | |  / \  / ___|| ____|
| |_) | |_| | / _ \ \___ \|  _|  
|  __/|  _  |/ ___ \ ___) | |___ 
|_|   |_| |_/_/   \_\____/|_____|
                                 
 ____ _____ ____  _   _  ____ _____ _   _ ____  _____ 
/ ___|_   _|  _ \| | | |/ ___|_   _| | | |  _ \| ____|
\___ \ | | | |_) | | | | |     | | | | | | |_) |  _|  
 ___) || | |  _ <| |_| | |___  | | | |_| |  _ <| |___ 
|____/ |_| |_| \_\\___/ \____| |_|  \___/|_| \_\_____|
                                                      
 _____ _   _ _   _  ____ _____ ___ ___  _   _ 
|  ___| | | | \ | |/ ___|_   _|_ _/ _ \| \ | |
| |_  | | | |  \| | |     | |  | | | | |  \| |
|  _| | |_| | |\  | |___  | |  | | |_| | |\  |
|_|    \___/|_| \_|\____| |_| |___\___/|_| \_|
                                              
*/

func computeFittingDphi_DMcircu_analytic(N, ud, r0tot, dactu)
/* DOCUMENT dphi = computeFittingDphi_DMcircu_analytic(N, ud, r0tot, dactu)

   Calcule analytiquement la fonction de structure de l'erreur de fitting
   sans passer par la TF du spectre, pour un DM a "maillage circulaire".
   Donne le mm resultat que
   W=Wfit; W(N/2+1,N/2+1) = 0.0;W(N/2+1,N/2+1) = -sum(W);Dphi = 2*abs(fft(W))*(uk^2);
   mais sans les erreurs de repliement, ni manque d'energie du spectre qui
   ne va pas jusqu'a +infty.
     
   SEE ALSO:
 */
{
    r = compute_rxy(N, ud);
    dphi = dphi_highpass(r,dactu) * r0tot^(-5./3);
    return dphi;
}

  
func compute_DphiBwFitting(r0tot, outerScale, V,wdir,Fe,tret,gain,BP,sgeom,mode=,verb=)
{
  
  // compute radial spatial frequencies in m^-1
  k = computeSpatialFreqRad(data.fourier.npix, data.fourier.uk, kx, ky);

  // mise a 0 correction DM : erreur de fitting
  // On definit la liste des indices <ndm> qui est l'ensemble des
  // frequences atteintes par le DM, et <northo> le complement
  //
  msk = defineDmFrequencyArea(k, kx, ky, sgeom, data.fourier.dactu );
  ndm    = where( msk );
  northo = where( !msk );
  
  
  // Wiener spectrum computation in m^-2
  Wiener = computeWienerSpectrum(k, data.fourier.npix, r0tot, outerScale);
  
  
  // Computation of Dphi of fitting error .............
  
  // We compute either Dphi (if DM geom is circular) or the phase
  // spectrum Wfit if DM has square arrangement
  Dfit = Wfit = [];
  if( sgeom=="circ" ) {
    Dfit = computeFittingDphi_DMcircu_analytic( data.fourier.npix, data.fourier.ud, r0tot, data.fourier.dactu);
    data.budget.fit = dphi2rms(Dfit);
  }
  else {
    Wfit = computeFittingWfit_DMsquare(Wiener, ndm, data.fourier.npix,verb=verb);
  }

  
  // Bandwidth error, computed in a different way. Here we use the
  // transfer function of the system, defined by a sampling frequency
  // and delay, assuming a pure close-loop integrator. The attenuation
  // of the transfer function Hcor(nu) is applied on all spatial
  // frequency k=nu/V (units equation is m^-1 = s^-1 / (m/s) ).
  // 
  Wbp = computeBpSpectrum(k,V,wdir,Fe,tret,gain,BP,Wiener,northo,mode=mode,verb=verb);
  
  
  // summing spectra ...........
  // including Wfit, if it's ever been computed
  W = Wbp;
  if( Wfit!=[] ){
    W += Wfit;
  }
  // normalisation pour valeur centrale spectre
  N = data.fourier.npix;
  W(N/2+1,N/2+1) = 0.0;
  W(N/2+1,N/2+1) = -sum(W);
  
  // Dphi computation in radians^2 at the good wavelength
  Dphi_bp = 2*abs(fft(W)) *  (data.fourier.uk^2);


  if( Dfit!=[] ) {
    Dphi = Dphi_bp + Dfit;
    return Dphi;
  } else {
    return Dphi_bp;
  }

  
}
/*
 ____  ____  _____ ____ _____ ____  _   _ __  __ 
/ ___||  _ \| ____/ ___|_   _|  _ \| | | |  \/  |
\___ \| |_) |  _|| |     | | | |_) | | | | |\/| |
 ___) |  __/| |__| |___  | | |  _ <| |_| | |  | |
|____/|_|   |_____\____| |_| |_| \_\\___/|_|  |_|
*/                                               



func computeSpatialFreqRad(N, uk, &kx, &ky)
{
  kx = ( indgen(N)-(N/2+1) ) * uk;   // on genere la composante X du vecteur freq spatiale k
  kx = kx(,-:1:N);                   // en 2 dims ...
  ky = transpose(kx);                // et pareil pour Y
  k = abs(kx,ky);
  k(N/2+1,N/2+1) = uk;           // valeur non-nulle pour ne pas calculer 0.00^(-11/3)=infini
  data.fourier.k = &k;
  data.fourier.kx = &kx;
  data.fourier.ky = &ky;
  return k;
}

func defineDmFrequencyArea(k, kx, ky, sgeom, dactu )
{
  fc = 0.5/dactu;
  
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


func computeFittingWfit_DMsquare(Wiener, ndm, N,verb=)
{
  
  Wfit = Wiener;
  Wfit(ndm) = 0.0;
  if(verb){
      write,format="Fit error from Wiener spectrum = %.4g nm\n", psd2rms(Wfit);
  }
  return Wfit;
}

func computeBpSpectrum(k,V,dir,Fe,tret,gain,BP,Wiener,northo,mode=,verb=)
/* DOCUMENT Wbp = computeBpSpectrum(k,V,dir,Fe,tret,gain,Wiener,northo)
     
   Bandwidth error. We use the transfer function of the system,
   defined by a sampling frequency and delay, assuming a pure
   close-loop integrator. The attenuation of the transfer function
   Hcor(nu) is applied on all spatial frequency k=nu/V (units equation
   is m^-1 = s^-1 / (m/s) ).

   SEE ALSO:
 */
{


  Wbp = 0*Wiener;
  cnh = data.learn.cnh;
  l0  = data.learn.l0h;
  vh = data.learn.vh;
  //Compensation of the model unaccuracy of the spatial covariance
  cnh *= data.turbu.r0ir^(-5/3.)/sum(cnh);
  L0eff = sum(l0^(5/3.)*cnh/sum(cnh))^(3/5.);
  l0  *= data.turbu.L0/sum(L0eff);

  for(l = 1;l<=data.learn.nl;l++){
    nu_l = k * vh(l)/sqrt(2);
    if(mode == "MOAO"){
      hcor = hcorMoao(nu_l, Fe, tret, gain, BP);
    }else if(mode == "SCAO"){
      hcor = hcorScao(nu_l, Fe, tret, gain, BP);
    }
    
    Wbp +=  hcor * 0.023 * cnh(l) * (k^2 + 1./l0(l)^2.)^(-11/6.);
  }

  
  N = data.fourier.npix;
  Wbp(N/2+1,N/2+1) = 0;
  Wbp(northo) = 0.00;
  
  if(verb){
    write,format="BW error from Wiener spectrum = %.4g nm\n", psd2rms(abs(Wbp));
  }

  return Wbp;
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



func compute_rxy(N, ud, &xx, &yy)
{
  xx = roll( ( indgen(N)-(N/2+1) ) * ud );
  xx = xx(,-:1:N);                
  yy = transpose(xx);             
  r = abs(xx,yy);         
  return r;
}




func psd2rms( W )
{
  fact = data.fourier.uk * data.camir.lambda_ir /2/pi * 1e9;
  return sqrt(sum(W)) * fact;
}

func dphi2rms(dphi)
{
  fact =  data.camir.lambda_ir /2/pi * 1e9;
  return sqrt(sum(abs(fft(-0.5*data.fourier.uld*dphi)/(numberof(dphi)-1.))))*fact;
}

func computeOthers(N, nmrms_others, uk, ndm, lambdaIR)
/* DOCUMENT 
     
   Autres sources d'erreur :
   cophasing  M1 ?????
   DM saturation, calibration and NCPA errors

   SEE ALSO:
 */
{
  Wothers = array(0.0, N, N);
  data.budget.ncpa = 0.00;

  fact = 2*pi* nmrms_others / uk / (lambdaIR*1e9);
  fact = fact^2;
  tot = dimsof(ndm)(0);          // number of pixels within dm space, avoid doing sum(...) to normalize
  Wothers(ndm) = fact / double(tot);
  data.budget.ncpa = psd2rms(Wothers);
  
  return Wothers;

}




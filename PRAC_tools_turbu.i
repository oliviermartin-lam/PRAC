func determineGlobalParameters(dataset,&p,&dp,arc=,verb=)
/* DOCUMENT noiseVarPix = determineGlobalParameters(dataset, &rzero, &lzero, &wspeed,ings=)
   
   Returns the variance of the noise in pix^2 for each slope of
   dataset (must be in pixel). Estimates rzero and lzero from
   projection of dataset on Zernike and the wind speed from the
   temporal autocorrelation.
     
 */
{
  // ............. Computation of noise ..............//
  noiseVar = getNoiseVar(dataset);

  // ............. Identification of OUTER SCALE AND R0 ......................//
  
  ings =  where(data.wfs.type==1);   // index of NGS WFS
  nngs = numberof(ings);
  nwfs  = data.nwfs;
  
  tmpR0L0 = du = array(0.0, nngs, 2);
  for(p=1; p<=nngs; p++) {
    icam        = ings(p);
    noise       = getNoiseZernike(noiseVar(slrange(icam)), icam, full=1, inputNoise=1,arc=arc);       // get noise in mic^2
    tmpR0L0(p,) = getr0L0( dataset(slrange(icam),), icam, noise,tmpdp, i0=i0,arc=arc);
    du(p,)      = tmpdp^2;
  }
  
  //..... WIND SPEED AND DIRECTION ......//
  v  = dir = dv = da = array(0.0, nwfs);
  w = array(0.0,nngs);
  delay = .1;k=1;
  for(p=1; p<=nwfs; p++) {
    if(data.wfs(p).type ==1){
      w(k)   = getWindSpeed(dataset(slrange(p),),p);
      k++;
    }
    tmp    = getWindFromCov(dataset(slrange(p),),delay,d1,d2,altlgs = data.wfs(p).lgsH);
    v(p)   = tmp(1);
    dir(p) = tmp(2);
    dv(p)  = d1^2;
    da(p)  = d2^2;
  }

  //derivating global turbulence parameters
  r0 = median(tmpR0L0(,1));
  L0 = median(tmpR0L0(,2));
  wspeed = avg(v);
  wdir = avg(dir);
  v2 = avg(w);
  //derivating uncertainties
  dr = sqrt(tmpR0L0(rms,1)^2/(nngs-1.) + du(avg,1));
  dl = sqrt(tmpR0L0(rms,2)^2/(nngs-1.) + du(avg,2));
  dws = sqrt(v(rms)^2/(nwfs-1.) + dv(avg));
  dwd = sqrt(dir(rms)^2/(nwfs-1) + da(avg));
  dv2 = w(rms)/sqrt(nngs-1.);
  //merging results
  p = [r0,L0,wspeed,wdir,v2];
  dp = [dr,dl,dws,dwd,dv2]

  return noiseVar;
}

/*

       ___                    _   _     ___  
 _ __ / _ \    __ _ _ __   __| | | |   / _ \ 
| '__| | | |  / _` | '_ \ / _` | | |  | | | |
| |  | |_| | | (_| | | | | (_| | | |__| |_| |
|_|   \___/   \__,_|_| |_|\__,_| |_____\___/ 
                                             
*/

func getr0L0( s, icam, noise,&du, i0=, i1=,arc=)
/* DOCUMENT getr0L0( s, icam, noise, i0=, i1=)
  
   <slopes> may be either
   - a slopes array, (72,1000..)
   - a vector of zernike variances
   
   <icam> is the camera number

   <noise> is either a scalar value (average noise on the WFS), or the
   noise on all zernike. In both cases it must be in mic^2. <noise>
   should be computed with the function noise_a_la_rico(s,icam,full=0/1).
   
   SEE ALSO: noise_a_la_rico(s, icam)
 */
{
  // Zernike reconstruction ...............
  if( dimsof(s)(1)==2 )
    varz = varZernike( s , icam ,arc=arc);     // in (nm rms)^2
  else
    varz = s;     // has to be in (nm rms)^2, just as the output of varianceZernike( s , icam );
  // conversion from (nm^2) to (rd rms)^2
  varz *=  (2*pi/(data.lambda_vis*1e9))^2;

  if(sum(varz<0)) return 0;
  // Noise propagation ....................
  noisez = 0.0;
  if( !is_void(noise) ) {
    if( isscalar(noise) ) {
      if( noise>0 ) {
        // noise, as computed by function noise_a_la_rico(s,icam), in mic^2.
        // <z_i^2> = $_j R_ij^2 <s_j^2>
        propj = ((*data.wfs(data.its).mrz)^2)(,sum);            // propag coeff on zernike
        noisez = noise * propj / propj(sum);     // noise variance on each zernike, mic^2
        noisez *= (2*pi/(data.lambda_vis*1e6))^2;                  // from mic^2 to radians^2 at 0.5 microns
      }
    } else {
      // noise, as computed by function noise_a_la_rico(s,icam,full=1), in mic^2.
      noisez = noise * (2*pi/(data.lambda_vis*1e6))^2;                  // from mic^2 to radians^2 at 0.5 microns
    }
    nn = where( varz>noisez );
    nn0 = where( varz<=noisez );
    
    if( is_array(nn) )
      varz(nn) -= noisez(nn);                // correction du bruit !!
    if( is_array(nn0) )
      varz(nn0) /= 100.;                     // division arbitraire par 100 quand bruit>signal
  }

  if(min(abs(varz))<1e-6) return 0;


  // i0 is the index where the fit will start
  if(is_void(i0))
    i0 = 3;     // to avoid tiptilt
  if(is_void(i1))
    i1 = 35;     // to avoid high orders
  
  // initial guess with D/r0=42 and D/L0=1
  Dr0L0 = data.tel.diam/[0.15, 20.0];
  res = lmfit(fitKolmo, [7,i0,i1], Dr0L0, log(varz(i0:i1)),stdev=1,tol=1e-7);
  du = (*res.stdev)*data.tel.diam/(Dr0L0)^2.;
  
  return data.tel.diam/abs(Dr0L0);
}


func varZernike( s , icam ,arc=)
/* DOCUMENT varzer = varZernike( s , icam )
     Compute variance of zernike polynomials.
     The variable <s> can either be a slope vector, or an array of slopes.
     For a single vector, the square of the zernike coeffs is returned.
     For a series, the centered mean-square value is returned.

     The variance array returned goes from TILT to Z35
   SEE ALSO:
 */
{ 
  // zernike reconstruction on zernike modes
  z = ((*data.wfs(icam).mrz)(,+)*s(+,)) ;
  // variance
  if(dimsof(s)(1)==1 || dimsof(s)(3)==1) {        // if s has 1 slope
    varz = z^2.;
  } else {
    varz = z(,rms)^2.;
  }
  units = (data.wfs(icam).unit)*1000.;
  if(arc) units = ((1/206265.)*(0.5*(data.tel.diam)*1e9));
  varz *= units^2;      // in nm^2
  return varz;
}


func zernikeSpectrum(radorder, D_r0, D_L0)
/* DOCUMENT zernikeSpectrum(radorder, D_r0, D_L0)

   <radorder> is the max radial order
   <D_r0> is the value of (D/r0),
   <D_L0> is the value of (D/L0).

   Returns the theoretical zernike spectrum of the turbulence, with
   D/r0 and D/L0 given as arguments, and seen+reconstructed though a
   7x7 subaperture SH.
   The influence of the 7x7 hartmann is simulated by multiplying the
   theoretical spectrum by coefficients (that have been calculated
   elsewhere) that simulate the aliasing impact.
   
   SEE ALSO:
 */
{
  varzer = zernikeVonKarmanSpectrum(((radorder+1)*(radorder+2))/2, abs(D_L0), abs(D_r0));

  zerfactor = [1.00127,1.00127,1.0089,1.00421,1.01433,1.05697,1.05697,1.0464,1.0464,1.07974,
               1.12548,1.03053,1.07697,1.05695,1.0647,1.0647,1.13497,1.13497,1.12401,1.12401,
               1.46191,1.32155,1.30285,0.752463,1.25776,0.911428,0.994586,2.69474,2.69474,
               2.18867,2.18867,1.3627,1.3627,1.39158,1.39158];

  varzer *= zerfactor;
   
  return varzer;
}

func zernikeVonKarmanSpectrum(jmax, dlo, dro)
/* DOCUMENT zernikeVonKarmanSpectrum(jmax, dlo, dro)
   
     <jmax> is the number of the last Zernike,
     <dlo> is D/L0
     <dro> is D/r0
     The function returns the theoretical spectrum of Zernike coeffs,
     with a given r0 and outer scale. It uses a formula given by Rodolph Conan in his phD thesis.

     The returned vector ranges from Z_2 (tilt) to Z_jmax.
     The piston variance is not returned (it could, but it's not ...).

     The function works for large L0>D/4, i.e. for 0<= (D/L0) < 4,
     it will diverge for L0<D/4.
     
   SEE ALSO:
 */
{
  local n,m;
  // Number of terms of the taylor expansion.
  // This number can be increased if higher accuracy is required,
  // e.g. for computing spectra for small L0<D/4, however the success
  // is not guaranteed.
  kmax=50;

  // mem allocation for result
  result=array(0.0, jmax-1);
  pidf0=pi*dlo; //produit pi*D*fo = pi*D/L0
  //  multiplicative factor  1.1628911592006026 :
  //  fact=2*gamma_R(11./6.);
  //  fact/=pi^(3./2.);
  //  fact*=((24./5)*gamma_R(6./5))^(5./6);
  dro53_116 = dro^(5./3) * 1.1628911592006026;
  // arrays n and m are outputs. They are the radial and azim order for polynoms 2 to jmax
  // n and m have (jmax-1) elements
  ordreNM,indgen(2:jmax),n,m;
  nmax = max(n);

  // The convergence of the series is not garanteed for values of D/L0<4 (i.e. pi*D/L0<4*pi)
  // Function may return something stupid when D/L0<4.
  // A test was previously setup so that the function returns 0.
  // Now (EG, 21 feb 2014) the function limits the L0 values (clipping) to L0>D/4
  // 
  if(pidf0>pi*4) {
    pidf0=pi*4;
  }
  // case D<L0
   for(in=1; in<=nmax; in++) {
    s = series_sum_diag(in,kmax,pidf0) * (in+1) * dro53_116;
    result(where(n==in)) = s;
   }
  return result;
    
}

func ordreNM(i, &n, &m)
/* DOCUMENT nm = ordreNM(i, n, m)
     Pour un zernike donne, defini par son indice i, calcule
     l'ordre azimutal m et l'ordre radial n.
 */
{
  n = int( (-1.+sqrt(8*(i-1)+1))/2.);
  p = (i-(n*(n+1))/2);
  k = n%2;
  m = int((p+k)/2)*2 - k;
}

func series_sum_diag(n,kmax,pidf0)
/* DOCUMENT series_sum_diag(n,kmax,pidf0)

   This function computes a sum of terms composed with Gamma
   functions, for the computation of zernike variances with finite
   L0. The formula is extracted from R. Conan thesis (equation 4.16
   page 116).
   The formula has been modified to use the relation
   Gamma(k+x)=(k-1+x).Gamma(k-1+x) and thus avoid to compute the Gamma
   functions at each iteration.
   
   SEE ALSO:
 */
{
  n2 = 2*n;

  // compute n! = Gamma(1+n)
  fn = 1.00;
  // compute (2+n1+n2)! using the previous results again
  for(i=2; i<=2+n2; i++) fn*=i;  // (2+n+n)! = Gamma(3+n+n)

  // computation of all gamma_R as a whole
  pregamma = gamma_R([1.5 + n, 5./6 - n, n-5./6, n+23./6]);

  // Initialisation
  uk = pregamma(1) * pregamma(2) / ((n+1) * fn);  // u0

  //  0.6494396447776356 = gamma_R(7/3.)*gamma_R(11./6)/gamma_R(17./6)
  vk = pregamma(3)*0.6494396447776356 / pregamma(4);
    
  pidf0_n253 = pidf0^(n2-5/3.);
  pidf0_2 = pidf0*pidf0;
  pidf0_2k = 1.00;   // preparing pidf0^(2*k), for k=0
  fk = 1.0;          // preparing k!, for k=0

  s = (uk * pidf0_n253 + vk);    // s(k=0)
  
  for(k=1; k<=kmax; k++) {
    fk *= k;    // k!
    pidf0_2k *= pidf0_2;   // computation of  pidf0^(2*k)
    
    uk *= ((0.5+n+k)*(k+n))/((2+n2+k)*(1+n+k)*(5./6-n-k));
    vk *= (k+4./3)*(k+5./6)/((n-5./6-k)*(n+17./6+k)*(11./6+k));
                                       
    tmps = (pidf0_n253 * uk + vk) * pidf0_2k / fk;
      
    if(k%2) s -= tmps;
    else  s += tmps;                                     
  }
  
  return s;
}


func fitKolmo(x,a)
{
  radorder = x(1);
  i0 = x(2);
  i1 = x(3);
  tmp = zernikeSpectrum(radorder, a(1), a(2))(i0:i1);
  return log(abs(tmp));
  // abs is not required, but allows to avoid bugs when running on crap data ..

}

/*
 _   _       _          
| \ | | ___ (_)___  ___ 
|  \| |/ _ \| / __|/ _ \
| |\  | (_) | \__ \  __/
|_| \_|\___/|_|___/\___|
                        
*/
func getNoiseZernike(s, icam, full=, inputNoise=,arc=)
/* DOCUMENT noise_a_la_rico(slopes, icam, full=, inputNoise=)
   
     Returns the average noise variance, in microns^2 (computed on the zernike).
     If <full>==1, returns the noise variance for each Zernike.

     The input data <s> may either be a slope array (72xN), or the
     array of the noise variance (in pix^2) for each slope.
     
     The noise variance is computed in pixels^2 on the slopes, and
     then propagated through the Zernike reconstruction matrix
     *rtc.mrz, and cumulated over the 36 Zernike.

     Returns microns^2.
     
   SEE ALSO:
 */
{
  extern rtc;

  if( inputNoise==1 )
    varNoise = s;  // assume that input data are noise variances, in pix^2
  else
    varNoise = getNoiseVar(s);  // returns the noise variance for EACH subaperture, in pixel^2
  
  // propagation of each individual subaperture noise
  units = data.wfs(icam).unit(-,);
  if(arc) units = ((1/206265.)*(0.5*(data.tel.diam)*1e6));
  zernoiseVar = (((*data.wfs(icam).mrz)^2)(,+) * varNoise(+)) * units^2;
  if( is_void(full) ) full=0;
  if( full==1 )
    return zernoiseVar;       // returns the full zernike variances, mic^2
  else
    return sum(zernoiseVar);  // returns the total propagated variance, mic^2
}


func getNoiseVar(s, zero=)
/* DOCUMENT noise_inPix2(s, zero=)
   
     Returns an array of the noise variance, in pixel^2,
     of each subaperture.
     Parabolic fit of the 3 1st points.
     
   SEE ALSO:
 */
{
  npt = dimsof(s)(3);           // number of frames in the slope set
  if( npt<5 )                   // meaningless result anyway
    return 0.0;
  sa = s(,avg);
  a0 = ((s(,1:npt)-sa)^2)(,avg);
  a1 = ((s(,1:npt-1)-sa)*(s(,2:npt)-sa))(,avg);
  a2 = ((s(,1:npt-2)-sa)*(s(,3:npt)-sa))(,avg);
  // parabolic fit
  varNoise = a0 - (4*a1-a2)/3.;
  if(zero) {
    nn = where(varNoise<0);
    if(is_array(nn)) varNoise(nn)=0.;
  }
  return varNoise;
}


/*
__        ___           _                         _ 
\ \      / (_)_ __   __| |___ _ __   ___  ___  __| |
 \ \ /\ / /| | '_ \ / _` / __| '_ \ / _ \/ _ \/ _` |
  \ V  V / | | | | | (_| \__ \ |_) |  __/  __/ (_| |
   \_/\_/  |_|_| |_|\__,_|___/ .__/ \___|\___|\__,_|
                             |_|                    
*/
func getWindSpeed(s, icam )
{

  npt = dimsof(s)(3);                                   // number of frames in the slope set
  if( npt<5 )     // meaningless result anyway
    return 0.0;
  // computes the autocorrelation
  s -= s(,avg);
  autocor = (fft( abs(fft(s,[0,1]))^2 ,[0,1]).re)(avg,);
  //  re-interpolate the first point : this allows us to avoid bias due to noise variance
  auto0 = autocor(2)*2.-autocor(3);          // linear interp
  auto0 = (autocor(2)*4.-autocor(3))/3.;     // parabolic interp
  if( auto0!=0 )
    autocor /= auto0;       // normalization
  else
    return 0.0;    // probably one has got this because all data==0
  npt/=2;
  k=2;
  while( autocor(k)>0.5 & k<npt) k++;   // search first point below 0.5
  // linear regression on 5 points around the half-max point
  if( (k-2)>3 & (k+2)<npt ) {
    y = autocor(k-2:k+2);
    x = indgen(k-2:k+2);
    tmp = (avg(x*x)-avg(x)^2.);
    a = (avg(y*x)-avg(x)*avg(y))/tmp;
    b = avg(y) - a*avg(x);
    if( a<0 )
      k = (0.5-b)/a;  // solve equation a*k+b = 0.5
    else
      k=0;
  }
  tau0 = k / data.rtc.Fe;   // tau0 in seconds

  // windspeed
  size_ssp = data.tel.diam/data.wfs(icam).sX;
  if( tau0>0 )
    windspeed = 1.61*size_ssp/tau0;// factor 1.61 has been "calibrated" by simulation
  else
    windspeed=0.0;
  return windspeed;
}

func getWindFromCov(s,delay,&dv,&da,altlgs=)
{
  if(is_void(delay)){delay = .1;}//temporal delta in s
  if(is_void(altlgs)) altlgs =0;
  s -= s(,avg);

  tau = int(delay*data.rtc.Fe)/data.rtc.Fe;
  if(tau == 0){
    tau   = 1/data.rtc.Fe;
    delay = 1/data.rtc.Fe;
    write,"\rDelay to short, set to 1 pixel shift";
  }
  //Computing the spatio-temporal covariance matrix
  Css = s(,+) * (roll(s,[0,1]*int(data.rtc.Fe*delay)))(,+)/dimsof(s)(0);

  //..... Fitting of the pupil shift .....//
  
  //managing intial guess
  data.learn.nl = 1.;
  data.learn.altitude = 0.;
  if(data.learn.cnh(1) == 0){data.learn.cnh(1) = 20;}
  if(data.learn.l0h(1) == 0){data.learn.l0h(1) = 10;}
  //managing index fit
  data.learn_fit.cnh(1) = 1;
  data.learn_fit.l0h(1) = 0;
  data.learn_fit.altitude(1) = 0;
  data.learn_fit.xshift(1) = 1;
  data.learn_fit.yshift(1) = 1;
  fitEstim = packcoeffs(data,indexFit);
  //fit the covariance
  F = calc_TTMatFilt_1(data.wfs(data.its).nssp);
  Css =  (F(,+)*Css(+,))(,+) * F(+,);
  takesDiag(Css) = 0;
  res = lmfit( covMat1WFS, data, fitEstim, Css, fit=indexFit, tol=1e-4,stdev=1);
  unpackcoeffs, fitEstim, cnh, alt,L0,tracking,xshift,yshift;

  //determining windspeed
  vx  = xshift(1)/tau;
  vy  = yshift(1)/tau;
  v   = abs(vx,vy);
  dir = 180*atan(vy/vx)/pi;

  //computing uncertainties
  sig = *res.stdev;
  sig = sig(where(sig!=0));
  dvx = sig(1)/tau;
  dvy = sig(2)/tau;
  dv  = abs(dvx,dvy);
  da  = 180/pi * (vx*dvy + vy*dvx)/v^2;
    
  return [v,dir]
}

func covMat1WFS(data,fitEstim)
/* DOCUMENT calcCaa_3(nbWfs, L0, r0, altitude)
     
   <nbWfs>      : scalar, number of wfss
   <L0>         : scalar, log10(outer scale in meters)
   <r0>         : array 1D, lenght = # of layers
   <altitude>       : array 1D, lenght = # of layers
   
   
   SEE ALSO:
 */
{
  altlgs = data.wfs(1).lgsH;

// get all the variables that were encapsulated inside fitEstim
  unpackcoeffs, fitEstim, cnh, alt, l0h, tracking, xshift, yshift;

  l0h = abs(l0h);
  if(l0h(1) >=100.) l0h(1)=100.;
  // reservation memoire
  covMat = array(0.,data.wfs(data.its).nssp, sum(data.wfs(data.its).nssp));

  //geometry
  dipl = compute_dipl(data.wfs(data.its).sX,data.tel.diam,0,0,0,0,0,0,obs=data.tel.obs);
 
  djql = compute_dipl(data.wfs(data.its).sX,data.tel.diam,0,0,0,0,xshift(1),yshift(1),obs=data.tel.obs);


  // subap size in the pupil plane
  sspSize = data.tel.diam/data.wfs(data.its).sX;
  
  covMat = Dij2covMat(dipl,djql,sspSize,sspSize,cnh(1),0,l0h(1),0,0,1,1);

  //removing the tip-tilt
  F = calc_TTMatFilt_1(data.wfs(data.its).nssp);
  covMat =  (F(,+)*covMat(+,))(,+) * F(+,);
  //nulifying the diagonal
  takesDiag(covMat) = 0;

  return covMat;
}


func getWindspeedProfile(&dv_h,&ddir_h,verb=)
/* DOCUMENT v =  getWindspeedProfile(dvh,ddirh,verb=1)

   Returns the retrieval of the turbulent profile performed by the L&A
   algorithm on the type- data_type (slopesdis, slopestl or datatomo)
   data set referenced by data_name (ex: "06h42m13s").
   
*/
{
  //keeping profile in memory
  learntmp = data.learn;
  nl     = data.learn.nl;
  cnhall = data.learn.cnh(1:nl);
  altall = data.learn.altitude(1:nl);
  l0hall = data.learn.l0h(1:nl);
  data.learn.diagonal= 0;
  data.learn.ttr = 0;
  data.learn_fit.tracking = 0;
  data.learn.tracking = [0,0,0];

  //..... Derivating the WFS measurements covariance matrix to be inverted .....//

  fitEstim = packcoeffs(data);
  Call = covMatModel(data, fitEstim,verb=verb);
  //adding noise and managing TT
  takesDiag(Call) += *data.turbu.varNoise;
  Call = handle_tilt_from_wfstype(Call);
  //determines number of modes to be filtered for the inversion
  condCaa = 300;
  mode = computeModesToBeFiltered(Call,data.its,condCaa = condCaa);
  //inversion of the full matrix
  Call_1 = inverse_mat(Call, mode);

  //managing the altitude covariance matrix computation
  nl = data.learn.nl;
  slopesdis = *data.slopes_dis;
  slopesdis -= slopesdis(,avg);
  
  delay = array(.15,data.learn.nl);
  delay(1) = .2;
  v_h = dir_h = dv_h = ddir_h = array(0.,nl);
  nwfs = data.nwfs;
  v = dir = dv = da = array(0.,nwfs);

  //Loop on layers
  for(l=1;l<=nl;l++){
    if(verb){
      write,format="Processing the %dth layer...\r",l;
    }
    // managing the data.learn struct
    data.learn.nl = 1;
    data.learn.cnh(1) = cnhall(l);
    data.learn.altitude(1) = altall(l);
    data.learn.l0h(1) = l0hall(l);
    
    //altitude covariance matrix computation;
    fitEstim = packcoeffs(data);
    Calt_hl = covMatModel(data,fitEstim,verb=verb,loworder=1);

    //Computing slopes from single layer
    R_hl = Calt_hl(,+)*Call_1(+,);
    slopes_hl = R_hl(,+)*slopesdis(+,);

    //Derivating the wind speed and directions
    for(p=1; p<=nwfs; p++){
      if(data.wfs(p).type == 1){
        tmp    = getWindFromCov(slopes_hl(slrange(p),),delay(l),d1,d2);
        v(p)   = tmp(1);
        dir(p) = tmp(2);
        dv(p)  = d1^2;
        da(p)  = d2^2;
      }
    }
    w   = where(v !=0);
    n   = numberof(w);
    v2   = v(w);
    dir2 = dir(w);
    dv2  = dv(w);
    da2  = da(w);
    //derivating global turbulence parameters
    v_h(l)    = avg(v2);
    dir_h(l)  = avg(dir2);
    //derivating uncertainties
    dv_h(l)   = sqrt(v2(rms)^2/(n-1.) + dv2(avg));
    ddir_h(l) = sqrt(dir2(rms)^2/(n-1.) + da2(avg));
  }
  
  //Getting back to the intial data struct
  data.learn = learntmp;
  data.learn.vh = v_h;
  data.learn.dirh = dir_h;
  data.uncertainties.vh = dv_h;
  data.uncertainties.dirh = ddir_h;
  
  return [v_h,dir_h];
}


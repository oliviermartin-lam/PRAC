func computePSDnoiseClosedLoop(Fe,tret,gain,BP,northo,mode,verb=)
/* DOCUMENT Wbp = computeBpSpectrum(k,V,dir,Fe,tret,gain,Wiener,northo)
     
   Bandwidth error. We use the transfer function of the system,
   defined by a sampling frequency and delay, assuming a pure
   close-loop integrator. The attenuation of the transfer function
   Hcor(nu) is applied on all spatial frequency k=nu/V (units equation
   is m^-1 = s^-1 / (m/s) ).

   SEE ALSO:
 */
{

  N = tel.nPixels;

  // Computing the PSD of off-axis noise propagated through reconstructor
  Cnn = *covMatrix.noise;
  if(mode == "MOAO"){
    sigNoise = computesNoisesError(takesDiag(Cnn),*rtc.R)(2);
  }else{
    s = (*rtc.slopes_dis)(slrange(rtc.its),);
    s -= s(,avg);
    sigNoise = sqrt(getNoiseZernike(s,rtc.its,arc=1))*1e3;
  }
  PSD_noise_ol  = (2*pi*sigNoise / tel.fourierPixSize / (cam.lambda*1e9))^2;
  
  nu = indgen(N)*Fe;
  
  if(mode == "MOAO"){
    hn = hboMoao(nu, Fe, tret, gain, BP);
  }else if(mode == "SCAO"){
    hn = hcorScao(nu, Fe, tret, gain, BP);
  }

  msk = array(1.,N,N);
  msk(northo) = 0;
  PSD_noise_cl =  sum(abs(hn)^2) * PSD_noise_ol * msk;

  PSD_noise_cl(N/2+1,N/2+1) = 0;  

  return PSD_noise_cl;
}

func getNoiseZernike(s, icam, full=, inputNoise=,arc=)
/* DOCUMENT noise_a_la_rico(slopes, icam, full=, inputNoise=)
   
     Returns the average noise variance, in microns^2 (computed on the zernike).
     If <full>==1, returns the noise variance for each Zernike.

     The input data <s> may either be a slope array (72xN), or the
     array of the noise variance (in pix^2) for each slope.
     
     The noise variance is computed in pixels^2 on the slopes, and
     then propagated through the Zernike reconstruction matrix
     *rtc.mrz, and cumulated over the 36 Zernike.

     Returns nm^2.
     
   SEE ALSO:
 */
{
  if( inputNoise==1 )
    varNoise = s;  // assume that input data are noise variances, in s^2 units
  else
    varNoise = getNoiseVar(s);  // returns the noise variance for EACH subaperture, in s^2 units
  
  // propagation of each individual subaperture noise
  units = wfs(icam).unit(-,);
  if(arc) units = .5*tel.diam*1e9/radian2arcsec;
  zernoiseVar = (((*sys.slopesToZernikeMatrix)^2)(,+) * varNoise(+)) * units^2;
  if( is_void(full) ) full=0;
  if( full==1 )
    return zernoiseVar;       // returns the full zernike variances, nm^2
  else
    return sum(zernoiseVar);  // returns the total propagated variance, nm^2
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

func calcMatCov_Noise( noise, spotWidth )
{
  Nslopes = data.Nslopes;
  Nsubap = Nslopes/2;

  diagnoise = array(0.0, 3, Nsubap);
  D = tomo.tel.diam;
  
  for(i=1; i<=data.nwfs; i++) {
    if( data.wfs(i).type && i!=data.its ) {
      // CASE OF NGS .....
      if( data.wfs(i).type==1 ) {
        rr = tsubrange(i);
        diagnoise(1:2, rr) = noise(i);
        diagnoise(3, rr) = 0.0;
      }
      
      // CASE OF LASER GUIDE STAR .....
      if( data.wfs(i).type==2 ) {
        nssp = data.wfs(i).sX;
        xy = compute_ri(nssp, D, data.tel.obs);
        xy = transpose(xy);  // now an array(double,2,N)
        
        // let's find where the laser launch telescope is .....
        llt = abs(data.wfs(i).x, data.wfs(i).y);
        if( llt==0 )
          llt = [0.,0.];   // LLT is at pupil center
        else
          llt = [data.wfs(i).x, data.wfs(i).y] / llt * D/2;  // LLT is at pupil border
        
        // coordinates difference between subapertures and LLT
        xy -= llt;

        // spot elongation .. in fact, it's just proportional to coord vector
        xy = 206265. * data.wfs(i).lgsdH*xy/data.wfs(i).lgsH^2.;   // extension at Fwhm, in arcsec
        lgsExt = abs( xy(2,), xy(1,) );   // lengh of the extension
        lgsTheta = atan( xy(2,), xy(1,) );   // angle of extension
        totalExt = abs(lgsExt, spotWidth); // lengh of the extension including seeing, laser size, ...

        // scaling the noise with square of extension. The noise
        // across the small axis is still noise(i).
        noiseLongAxis = noise(i) * (totalExt / spotWidth)^2;
        
        s2x = noiseLongAxis * cos(lgsTheta)^2 + noise(i)*sin(lgsTheta)^2;
        s2y = noiseLongAxis * sin(lgsTheta)^2 + noise(i)*cos(lgsTheta)^2;
        sxy = (noiseLongAxis-noise(i)) * sin(lgsTheta)*cos(lgsTheta);

        // putting the result at the right place
        rr = tsubrange(i);
        diagnoise(1,rr) = s2x;
        diagnoise(2,rr) = s2y;
        diagnoise(3,rr) = sxy;
      }
    }
  }

  return diagnoise;
}







func debugcalcMatCov_Noise( i, diagnoise, width= )
{
  
  Nslopes = sum(tomo.wfs.Nslopes);
  Nsubap = sum(tomo.wfs.Nsubap);

  D = tomo.tel.diam;

  rr = tsubrange(i);
  diagos = diagnoise(,rr);
  
  nssp = tomo.wfs(i).nssp;
  xy = vectmesFab(nssp, D, tomo.tel.obs, 0 , 0, nssp, 0);

  n = dimsof(xy)(2); // number of subapertures

  for(i=1;i<=n;i++) {
    plot_ellipse, diagos(1,i), diagos(2,i), diagos(3,i), xy(i,1), xy(i,2), width=width ;
  }
        
}


func plot_ellipse( sx2, sy2, sxy, x0, y0, width= )
{
  ds = sx2 - sy2;
  if( ds==0 ) th = -pi/4.;
  else th = atan(2*sxy/ds)/2.;
    
  if( th==0 )
    sms=0.;
  else
    sms = 2*sxy/sin(2*th);
  sps = sx2 + sy2;

  SX = (sps+sms)/2;
  SY = (sps-sms)/2;
  t=span(0,2*pi,40);
  x = SX*cos(t);
  y = SY*sin(t);
  xx = x*cos(th)-y*sin(th);
  yy = x*sin(th)+y*cos(th);
  plg, yy+y0, xx+x0, marks=0, color="blue", width=width; 
}





func add_diagonal( mat, noisediag )
/* DOCUMENT add_diagonal, mat, noisediag
     
   SEE ALSO:
 */
{
  for(k=1; k<=tomo.naso; k++) {
    N = tomo.wfs(k).Nsubap;   // nber of subapertures
    rr = tsubrange(k);        // range where data are stored in <noisediag> for this wfs
    i0 = tslindex(k);         // 1st index where slopes are stored for this wfs

    // add the x-variance on matrix bocks of slopes in x
    wfsnoise = noisediag(1,rr);
    for(i=1; i<=N; i++) mat(i+i0-1,i+i0-1) += wfsnoise(i);
    
    // idem in y
    wfsnoise = noisediag(2,rr);
    for(i=1; i<=N; i++) mat(i+i0-1+N,i+i0-1+N) += wfsnoise(i);

    // and 2 diagonals for covariances
    wfsnoise = noisediag(3,rr);
    for(i=1; i<=N; i++) {
      mat(i+i0-1+N,i+i0-1) += wfsnoise(i);
      mat(i+i0-1,i+i0-1+N) += wfsnoise(i);
    }
  }
}

func noiseTiltLGS( noise, mat )
{
  n = dimsof(mat)(0);
  for(i=1; i<=tomo.naso; i++) {
    if( tomo.wfs(i).type==2 && i!=tomo.its ) {
      rr = tslindex(i):tslindex(i)+tomo.wfs(i).Nsubap-1;
      mat(rr,rr) += noise;
      rr = tslindex(i)+tomo.wfs(i).Nsubap:tslindex(i)+tomo.wfs(i).Nslopes-1;
      mat(rr,rr) += noise;
    }
  }
}

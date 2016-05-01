func determineGlobalParameters(dataset,&p,&dp,arc=,verb=)
/* DOCUMENT noiseVarPix = determineGlobalParameters(dataset, &rzero, &lzero, &wspeed,ings=)
   
   Returns the variance of the noise in arcsec^2 for each slope of
   dataset (must be in pixel). Estimates rzero and lzero from
   projection of dataset on Zernike and the wind speed from the
   temporal autocorrelation.
     
 */
{
  // ............. Computation of noise ..............//
  noiseVar = getNoiseVar(dataset);

  // ............. Identification of OUTER SCALE AND R0 ......................//
  
  ings =  where(wfs.type==1);   // index of NGS WFS
  nngs = numberof(ings);
  nwfs  = rtc.nWfs;
  
  tmpR0L0 = du = array(0.0, nngs, 2);
  v = array(0.0,nngs);
  for(p=1; p<=nngs; p++) {
    icam        = ings(p);
    // get noise in nm^2
    noise       = getNoiseZernike(noiseVar(slrange(icam)), icam, full=1, inputNoise=1,arc=arc);
    // get r0 and L0 from Zernike variance fitting
    tmpR0L0(p,) = getr0L0( dataset(slrange(icam),), icam, noise,tmpdp, i0=3,arc=arc);
    du(p,)      = tmpdp^2;
    //get wind velocity
    v(p)        = getWindSpeed(dataset(slrange(icam),),icam,noiseVar(slrange(icam)));
  }

  //derivating global turbulence parameters
  r0 = median(tmpR0L0(,1));
  L0 = median(tmpR0L0(,2));
  wspeed = avg(v);
  //derivating uncertainties
  if(nngs !=1){
    dr = sqrt(tmpR0L0(rms,1)^2/(nngs-1.) + du(avg,1));
    dl = sqrt(tmpR0L0(rms,2)^2/(nngs-1.) + du(avg,2));
    dv = v(rms)/sqrt(nngs-1.);
  }else{
    dr = sqrt(tmpR0L0(rms,1)^2 + du(avg,1));
    dl = sqrt(tmpR0L0(rms,2)^2 + du(avg,2));
    dv = v(rms);
  }
  //merging results
  p = [r0,L0,wspeed];
  dp = [dr,dl,dv]

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
  varz *=  (2*pi/(atm.lambda*1e9))^2;

  if(sum(varz<0)) return 0;
  // Noise propagation ....................
  noisez = 0.0;
  if( !is_void(noise) ) {
    if( isscalar(noise) ) {
      if( noise>0 ) {
        // noise, as computed by function noise_a_la_rico(s,icam), in mic^2.
        // <z_i^2> = $_j R_ij^2 <s_j^2>
        propj = ((*sys.slopesToZernikeMatrix)^2)(,sum); // propag coeff on zernike
        noisez = noise * propj / propj(sum);     // noise variance on each zernike, nm^2
        noisez *= (2*pi/(atm.lambda*1e9))^2;           // from nm^2 to radians^2 at 0.5 microns
      }
    } else {
      // noise, as computed by function noise_a_la_rico(s,icam,full=1), in nm^2.
      noisez = noise * (2*pi/(atm.lambda*1e9))^2;    // from nm^2 to radians^2 at 0.5 microns
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
  Dr0L0 = tel.diam/[0.15, 20.0];
  res = lmfit(fitKolmo, [7,i0,i1], Dr0L0, log(varz(i0:i1)),stdev=1,tol=1e-7);
  du = (*res.stdev)*tel.diam/(Dr0L0)^2.;
  
  return tel.diam/abs(Dr0L0);
}


/*
__        ___           _                         _ 
\ \      / (_)_ __   __| |___ _ __   ___  ___  __| |
 \ \ /\ / /| | '_ \ / _` / __| '_ \ / _ \/ _ \/ _` |
  \ V  V / | | | | | (_| \__ \ |_) |  __/  __/ (_| |
   \_/\_/  |_|_| |_|\__,_|___/ .__/ \___|\___|\__,_|
                             |_|                    
*/
func getWindspeedProfile(&dv_h,&ddir_h,verb=)
/* DOCUMENT v =  getWindspeedProfile(dvh,ddirh,verb=1)

   Returns the retrieval of the turbulent profile performed by the L&A
   algorithm on the type- data_type (slopesdis, slopestl or datatomo)
   data set referenced by data_name (ex: "06h42m13s").
   
*/
{
  //keeping profile in memory
  learntmp = learn;
  nl       = learn.nl;
  cnhall   = learn.cnh(1:nl);
  altall   = learn.altitude(1:nl);
  l0hall   = learn.l0h(1:nl);
  
  learn.diagonal = 0;
  learn.ttr      = 0;
  learn.tracking = 0;


  //..... Derivating the WFS measurements covariance matrix to be inverted .....//

  fitEstim = packcoeffs(learn);
  Call = covMatModel(learn, fitEstim,verb=verb);
  //adding noise and managing TT
  Call += *covMatrix.noise;

  //determines number of modes to be filtered for the inversion
  condCaa = 30.;
  mode = computeModesToBeFiltered(Call,rtc.its,condCaa = condCaa);
  //inversion of the full matrix
  Call_1 = inverse_mat(Call, mode);

  //managing the altitude covariance matrix computation
  nl = learn.nl;
  slopesdis = *rtc.slopes_dis;
  slopesdis -= slopesdis(,avg);

  //initialiazing vectors
  v_h  = dv_h  = array(0.,nl);
  ings =  where(wfs.type==1);   // index of NGS WFS
  nngs = numberof(ings);

  //Loop on layers
  for(l=1;l<=nl;l++){
    if(verb){
      write,format="Processing the %dth layer...\r",l;
    }
    // managing the data.learn struct
    learn.nl = 1;
    learn.cnh(1) = cnhall(l);
    learn.altitude(1) = altall(l);
    learn.l0h(1) = l0hall(l);
    
    //altitude covariance matrix computation;
    if(learn.altitude(1) >50)
      learn.tracking = [0,0,0];
    
    fitEstim = packcoeffs(learn);
    Calt_hl = covMatModel(learn,fitEstim,verb=verb);

    //Computing slopes from single layer
    R_hl = Calt_hl(,+)*Call_1(,+);
    slopes_hl = R_hl(,+)*slopesdis(+,);
    noiseVar = getNoiseVar(slopes_hl);

    //Derivating the wind speed
    v = array(0.0,nngs);
    for(p=1; p<=nngs; p++) {
      icam   = ings(p);
      v(p)   = getWindSpeed(slopes_hl(slrange(icam),),icam,noiseVar(slrange(icam)));
    }
    v_h(l)  = avg(v);
    if(nngs!=1)
      dv_h(l) = v(rms)/sqrt(nngs-1.);
    else
      dv_h(l) = v(rms)
  }
  
  //Getting back to the intial data struct
  learn = learntmp;
  atm.vh(1:nl)  = v_h;
  
  return vh_profile;
}

func getWindSpeed(s, icam,varNoise,&dir )
{

  // computes the autocorrelation
  s      -= s(,avg);
  nf     = dimsof(s)(3);    
  autocor = autocorrelation(s);
  autocor = fft(abs(fft(s,[0,1]))^2,[0,1]).re/nf^2;
    
  //denoising
  autocor(,1) -= varNoise;
  autocor     /= autocor(,max);
  a_x = (autocor(1:wfs(icam).nValidSubAp,))(avg,);
  a_y = (autocor(wfs(icam).nValidSubAp+1:,))(avg,);
  
  tau0_x = getFWHM(a_x,rtc.Fe);
  tau0_y = getFWHM(a_y,rtc.Fe);

  // windspeed
  size_ssp = tel.diam/wfs(icam).nLenslet;

  vx = 0.791*size_ssp/tau0_x;// factor 0.791 has been calibrated by yao simulation
  vy = 0.791*size_ssp/tau0_y;

  dir = atan(vy/vx)*180/pi;
  
  return abs(vx,vy);
}









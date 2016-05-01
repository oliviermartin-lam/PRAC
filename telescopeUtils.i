func pupilGeometry(N)
/* DOCUMENT [x,y] = pupilGeometry(N)

   Returns NxN array of x and y axis normalized
   for varying from -1 to 1.
*/
{
  x = (indgen(N) - (N+1.)/2.)(,-:1:N)/double(N/2.);
  //x = (indgen(N) - N/2.)(,-:1:N)/double(N/2.);
  y = transpose(x);
  return [x,y];
}
func pupilFunction(obs,N)
/* DOCUMENT P = pupilFunction(obs,pixSize)

   Returns a circular pupil function for a  telescope
   with a central obstruction obs, over N pixels.

    SEE ALSO: pupilGeometry
 */
{
  //defining geometry
  tmp = pupilGeometry(N);
  x   = tmp(,,1);
  y   = tmp(,,2);
  r   = abs(x,y);
  return r<=1.0 & r>obs;
}
func enlargeSupport(P,nTimes)
/* DOCUMENT P2 = enlargeSupport(P,nTimes)

   Enlarges by nTimes the size of function P
   in adding zero.
 */
{
  n  = dimsof(P)(0);
  P2 = array(0.,n*nTimes,n*nTimes);
  ri = n*(nTimes-1)/2+1;
  rf = n*(nTimes+1)/2;
  P2(ri:rf,ri:rf) = P;
  return P2;
}
func interpolateOTF(OTF,N)
/* DOCUMENT OTF2 = interpolateOTF(OTF,N)

   Interpolates the OTF function over
   a N x N support.
 */
{
  //Initial geometry
  N1  = dimsof(OTF)(0);
  tmp = pupilGeometry(N1);
  x1  = tmp(,,1);
  y1  = tmp(,,2);

  //Wished geometry
  tmp = pupilGeometry(N);
  x2  = tmp(,,1);
  y2  = tmp(,,2);

  //Interpolation
  return interp2(y2,x2,OTF,y1,x1);

}
func pracTelescopeOTF(obs,N,norm=)
/* DOCUMENT OTF = pracTelescopeOTF(obs,N)

  Returns the telescope OTF from the pupil autocorrelation function.
  The OTF is normalized to get a psf with PSF(0)=1.00 when diffraction-limited.
     
   SEE ALSO: pupilFunction, enlargeSupport,interpolate
 */
{
  //Computation of pupil function
  P       =  pupilFunction(obs,N);
  //Enlarging the pupil size
  P2      =  enlargeSupport(P,tel.nTimes);
  //Computing the OTF
  telOTF  = fft(abs(fft(P2))^2).re / numberof(P2)^2;
  // normalization
  if(norm)
    telOTF /= sum(telOTF);
  
  return telOTF;
}

func polarFromCartesian(x,y)
/* DOCUMENT [rho,theta] = polarFromCartesian(x,y)

   Transforms a Cartesian coordinate system into a polar one.
 */
{
  rho   = abs(x,y);
  theta = 0*rho;

  if(numberof(x) == 1){
    theta = atan(y/x);
    if(x>0 & y<0){
      theta += 2*pi;
    }else if(x<0){
      theta += pi;
    }else if(x==0 & y>0){
      theta = pi/2.;
    }else if(x==0 & y<0){
      theta = 3*pi/2.;
    }
  }else{
    w0        = where(x>0);
    if(is_array(w0))
      theta(*)(w0) = atan(y(*)(w0)/x(*)(w0));
    w0        = where(x>0 & y<0);
    if(is_array(w0))
      theta(*)(w0) = atan(y(*)(w0)/x(*)(w0)) + 2*pi;
    w0        = where(x<0);
    if(is_array(w0))
      theta(*)(w0) = atan(y(*)(w0)/x(*)(w0)) + pi;
    w0        = where(x==0 & y>0);
    if(is_array(w0))
      theta(*)(w0) = pi/2.;
    w0        = where(x==0 & y<0);
    if(is_array(w0))
      theta(*)(w0) = 3*pi/2.;
  }
  return [rho,theta];
}

func circularAveragePsf(im)
/* DOCUMENT imavg = circularAveragePsf(im)

  Returns the radial average of im.
*/
{
  // Geometry
  N   = dimsof(im)(0);
  x   = (indgen(N) -N/2-1)(,-:1:N);
  y   = transpose(x);
  tmp = polarFromCartesian(x,y);
  rho = tmp(,,1);
  t   = tmp(,,2); 
  rmax = max(rho);
  
  // Initialisation
  nt = 32;nr = N/2;//nr = 128;
  imAvg = npx = array(0.,nr,nt);
 
  for(i=1;i<=nt;i++){
    if(i == 1){
      msk_t = t < 2*pi/double(nt);
    }else{
      msk_t = (t >= (i-1)*2*pi/double(nt) ) & (t < i*2*pi/double(nt));
    }
    for(k=1;k<=nr;k++){
      if(k == 1){
        msk_r = rho < rmax/double(nr);
      }else{
        msk_r = rho >= (k-1)*rmax/double(nr) & rho < k*rmax/double(nr);
      }

      mskIm = im * msk_t * msk_r;
      w     = where(mskIm !=0);
      if(is_array(w)){
        imAvg(k,i) = avg(mskIm(w));
        npx(k,i)   = 1.;
      }
    }
  }

  w     = where(npx(,sum) !=0);
  imAvg = imAvg(w,sum)/npx(w,sum);
  
  return imAvg;
}



func staticOTF(&statPSF,&statSR,mode=,deconvolve=)
/* DOCUMENT OTF = staticOTF(statPSF,statSR)

   Calculates the static OTF by the e^(i x phi)
   auto-correlation, where phi is the static phase
   determined by the projection of the TS measurements
   into the mirror basis.
*/
{

  if(is_void(deconvolve))
    deconvolve = 1;

  if(is_void(mode))
    mode = "dm";

  //defining geometry
  N   = tel.nPixels;
  tmp = pupilGeometry(N);
  x   = tmp(,,1);
  y   = tmp(,,2);
    
  //computes static aberrations in arcsec
  res    = (*rtc.slopes_res)(slrange(rtc.its),);
  stats  = res(,avg);//arcsec
      
  // Derivating the phase
  if(mode == "dm"){

    mc  = *rtc.mcsky; // V/arcsec
    V   = mc(,+) * stats(+); //V
    Mi  = *dm.modes;
    phi = Mi(,,+) * V(+) * (2*pi/cam.lambda);
    
  }else if(mode == "zernike"){

    a      = (*sys.slopesToZernikeMatrix)(,+) * stats(+);
    a      = grow(0.,a); // adding piston;
    a     *=  pi*tel.diam/(radian2arcsec*cam.lambda);
    nmodes = numberof(a);
  
    tmp   = polarFromCartesian(x,y);
    rho   = tmp(,,1);
    theta = tmp(,,2);
  
    //loop on Zernike modes for CPA
    phi  = array(0.,N,N);
    for(i=2;i<=nmodes;i++){
      phi += a(i) *  computeZernikePolynomials(rho,theta,i);
    }
  }
  
  //Pupil function and Electric field
  P        = pupilFunction(tel.obs,N);
  E        = P*exp(1i*phi);
  //Doubling the size of the support
  P2       = enlargeSupport(P,tel.nTimes);
  realE2   = enlargeSupport(E.re,tel.nTimes);
  imE2     = enlargeSupport(E.im,tel.nTimes);
  //Autocorrelation
  autoP    = autocorrelation(P2);
  statOTF  = autocorrelation(realE2+1i*imE2);

  //computing statistics
  statSR = sum(statOTF)/sum(autoP);
  
  if(deconvolve){
    statOTF /= (autoP + 1e-15);
    statOTF *= autoP/max(autoP)>1e-5;
  }
  
  //normalization
  statPSF  = fft(statOTF).re;
  
  //Interpolating to the good resolution
  //statOTF = interpolateOTF(statOTF,tel.nPixels);

  return statOTF;
  
}

func changeResolution(otf,initSize,finalSize,finalFov)
/* DOCUMENT OTF = changeResolution(otf,initSize,finalSize,finalFov)

   Makes a zero padding or interpolation on the OTF to get
   the wished angular resolution in the direct space.
*/
{

  if(initSize == finalSize)
    return otf;

  // Setting the resolution
  N1   = dimsof(otf)(0);
  N2   = int(N1*initSize/finalSize);
  otf2 = array(0.,N2,N2);
  
  if(N2 > N1){
    xi = (N2 - N1)/2 + 1;
    xf = (N2 + N1)/2;
    otf2(xi:xf,xi:xf) = otf;
  }else{
    xi   = (N1 - N2)/2 + 1;
    xf   = (N2 + N1)/2;
    otf2 = otf(xi:xf,xi:xf);
  }

  if(finalFov){
    //Interoplating the OTF to get the wished fov
    N3   = int(finalFov/finalSize);
    otf2 = interpolateOTF(otf2,N3); 
  }
  
  return otf2;
}









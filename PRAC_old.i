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

func getWindFromCov(s,delay,&dv,&da,verb=)
{
  if(is_void(delay)){delay = .1;}//temporal delta in s
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
  tmptrack = data.learn.tracking;
  data.learn.tracking = 0;
  //managing index fit
  data.learn_fit.cnh      = 0;
  data.learn_fit.l0h      = 0;
  data.learn_fit.altitude = 0;
  data.learn_fit.xshift(*data.rtc.ptrlistOffAxisNgs) = 1;
  data.learn_fit.yshift(*data.rtc.ptrlistOffAxisNgs) = 1;
  data.learn.diagonal = 1;
  fitEstim = packcoeffs(data,indexFit);
  //fit the covariance
  F = calc_TTMatFilt_1(data.Nslopes);
  Css =  (F(,+)*Css(+,))(,+) * F(+,);
  takesDiag(Css) = 0;
  
  res = lmfit_Learn( covMatModel, data, fitEstim, Css, fit=indexFit, tol=1e-4,stdev=1,verb=verb);
  unpackcoeffs, fitEstim, cnh, alt,L0,tracking,xshift,yshift;

  //determining windspeed
  ws = where(xshift!=0);
  xs = xshift(ws);
  ys = yshift(ws);
  vx  = xs/tau;
  vy  = ys/tau;
  v   = abs(vx,vy);
  dir = 180*atan(vy/vx)/pi;

  //computing uncertainties
  sig  = *res.stdev;
  w    = where(sig!=0);
  sig  = sig(w);
  n    = numberof(w);
  sigx = sig(1:n/2);
  sigy = sig(n/2+1:);
  dvx  = sigx/tau;
  dvy  = sigy/tau;
  dv   = abs(dvx,dvy);
  da   = 180/pi * (vx*dvy + vy*dvx)/v^2;

  data.learn.tracking = tmptrack;
  return [v,dir]
}

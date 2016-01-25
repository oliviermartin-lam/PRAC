/*
 _____                       __  __       _        _               
| ____|_ __ _ __ ___  _ __  |  \/  | __ _| |_ _ __(_) ___ ___  ___ 
|  _| | '__| '__/ _ \| '__| | |\/| |/ _` | __| '__| |/ __/ _ \/ __|
| |___| |  | | | (_) | |    | |  | | (_| | |_| |  | | (_|  __/\__ \
|_____|_|  |_|  \___/|_|    |_|  |_|\__,_|\__|_|  |_|\___\___||___/                                  
 */

func computesCeeMatrix(obsmode,verb=)
{

  //retrieves the noise covariance matrix in arcsec^2
  Cnn = *covMatrix.noise;
  
  if(obsmode == "MOAO"){
    //turbulence contribution. We now have to include, noise, aliasing and tracking contributions from off-axis WFS through rtc.R into covOL
    covOL = *covMatrix.parallel;
    indOff = norange(rtc.its);
    
    //adding aliasing contribution from off-axis WFS
    covOL(indOff,indOff) += (*covMatrix.aliasing)(indOff,indOff);
    //adding the time-filtered off-axis noise contribution
    fact = filteringNoiseFactor(rtc.loopGain,rtc.delay*rtc.Fe,rtc.obsMode);
    covOL(indOff,indOff) += fact*Cnn(indOff,indOff);
    //adding tracking
    covOL(indOff,indOff) += (*covMatrix.tracking)(indOff,indOff);

    //derinving residue tomographic error
    Cee = computeCovSlopesError(covOL,*rtc.R);
    

  }else if(obsmode == "SCAO"){
    
    if(verb){
      write,"\rComputing the aliasing and noise injected into the loop..."
    }
    //noise propgated through the loop
    fgt  = filteringNoiseFactor(rtc.loopGain,rtc.delay*rtc.Fe,rtc.obsMode);
    Cee  =  fgt*Cnn(slrange(rtc.its),slrange(rtc.its));
    Cee += (*covMatrix.aliasing)(slrange(rtc.its),slrange(rtc.its));
  }
  
  return Cee;
}


 
func computeCovSlopesError(mat,R)
/* DOCUMENT eps = computeTomoError_cpu(mat,R)
   
     Computes the cov matrix of the residual error.
     
     Cee = computeTomoError_cpu(mat,R);

     The input covariance matrix <mat> is in arcsec².
     The output covariance matrix <Cee> is in arcsec², to be used
    
     SEE ALSO:
 */
{
  
  rr = 1:tslindex(1-rtc.nWfs);  // WFS 1 to last before TS
  tsr = tslrange(rtc.its);
  cpp = mat(tsr,tsr);
  cpm = mat(tsr,rr);
  cmm = mat(rr,rr);
  mat = [];
  
  // Computation of the cov matrix of residual phase error
  tmp = R(,+) * cpm(,+);
  eps = cpp - tmp - transpose(tmp) + R(,+) * (cmm(,+)*R(,+))(+,);
  
  return eps;
}

func computesTomoError(MatCov, R, TScam,full=)
/* DOCUMENT sigmatomo = computeTomoError(MatCov, R, TScam, verbose=,full=)

  Computes the tomographic error from the fitted covariance matrix
  expressed in arcsec^2, and from the reconstructor. This tomographic
  error is the minimal MMSE error plus the truth sensor's noise and
  the noise propagated through the reconstructor R.

  Olivier Martin
*/
{
  if(!is_void(R)){
    if(is_void(TScam)) TScam = rtc.its;
    YY = MatCov(slrange(TScam),slrange(TScam));
    XY = MatCov(norange(TScam),slrange(TScam));
    XX = MatCov(norange(TScam),norange(TScam));
          
    //Calcul erreur minimale avec la covariance du TS
    Cee = YY - XY(+,)*R(,+) - R(,+)*XY(+,) + R(,+)*(XX(,+)*R(,+))(+,);
  }else{
    Cee = MatCov(slrange(TScam),slrange(TScam));
  }
  
  //Derivation on the Zernike
  mrz = *sys.slopesToZernikeMatrix;
  Czz = mrz(,+)*(Cee(,+)*mrz(,+))(+,);
  nzz = dimsof(Czz)(0);
  // units ...
  units = 0.5*tel.diam*1e9/radian2arcsec;
  diagonale = takesDiag(Czz)*units^2;
  
  if(full){
    return SQRT(diagonale);
  }else{
    return SQRT((sum(diagonale)));
  }
}


/*
  ____                           _              
 / ___| ___  ___  _ __ ___   ___| |_ _ __ _   _ 
| |  _ / _ \/ _ \| '_ ` _ \ / _ \ __| '__| | | |
| |_| |  __/ (_) | | | | | |  __/ |_| |  | |_| |
 \____|\___|\___/|_| |_| |_|\___|\__|_|   \__, |
                                          |___/ 
*/

func getValidSubapArray( nLenslet, rext, rint )
/* DOCUMENT valid = getValidSubapArray( nssp, rext, rint )

   Example :
   > getValidSubapArray( 5, 1, 0.285 )
   [0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,0]
     
   SEE ALSO:
 */
{
  // tip-tilt sensor only.
  if(nssp==1){
    return [1];
  }
  // to avoid some bug that eliminates useful central subapertures when obs=0.286
  if(nssp==7 && (rint>0.285 && rint<0.29)){
    rint=0.285;
  }
  
  x=span(-1,1,nLenslet+1)(zcen)(,-:1:nLenslet);
  y=transpose(x);
  r=sqrt(x*x+y*y);
  valid2d = (r<rext & r>=rint);
  valid = valid2d(*);
  return valid;
}

func computeRi(nLenslet,diamtel,obs,xshift,yshift,magnification,theta)
/* DOCUMENT vectmesFab(nLenslet, diamtel, obs)
   
     
   SEE ALSO:
 */
{
  rtel = diamtel/2.;
  x = span(-rtel, rtel, nLenslet+1)(zcen)(,-:1:nLenslet);
  y = transpose(x);
  nn = where(getValidSubapArray( nLenslet, 1.0, obs ));
  // taking magnification factor into account
  x *= magnification;
  y *= magnification;
  // taking rotation into account
  th = theta * pi/180.;
  x_ = x*cos(th) - y*sin(th);
  y_ = x*sin(th) + y*cos(th);
  // taking pupil offsets (misregistration deviations) into account
  x_ += xshift;
  y_ += yshift;

  return [x_(*)(nn),y_(*)(nn)];

}


func computeDipl(nLenslet, D, h, x, y, Hlgs, xshift,yshift,magnification,theta,obs=)
/* DOCUMENT mefNGSorLGS(np, D, h, x, y, XPup, YPup, diamPup, thetaML, Hlgs, obs=)
     
   SEE ALSO:
 */
{
  if(is_void(obs)){
    obs=0.;
  }
  n = numberof(x);
  if(numberof(y)!=n){
    exit,"problem !! !";
  }
  
  dipl = computeRi(nLenslet,D,obs,xshift,yshift,magnification,theta)(..,-:1:n);
  
  if( Hlgs ) {  // if non-zero, the guide star is a LGS
      dipl(,1,) = dipl(,1,)*(Hlgs-h)/Hlgs + h*x;
      dipl(,2,) = dipl(,2,)*(Hlgs-h)/Hlgs + h*y;
  }else{
      dipl(,1,) += h*x;
      dipl(,2,) += h*y;
  }
  return dipl;
}


/*
  ____                     _                      
 / ___|_____   ____ _ _ __(_) __ _ _ __   ___ ___ 
| |   / _ \ \ / / _` | '__| |/ _` | '_ \ / __/ _ \
| |__| (_) \ V / (_| | |  | | (_| | | | | (_|  __/
 \____\___/ \_/ \__,_|_|  |_|\__,_|_| |_|\___\___|

 __  __       _        _               
|  \/  | __ _| |_ _ __(_) ___ ___  ___ 
| |\/| |/ _` | __| '__| |/ __/ _ \/ __|
| |  | | (_| | |_| |  | | (_|  __/\__ \
|_|  |_|\__,_|\__|_|  |_|\___\___||___/
                                       
*/

func dij2covMat(dipl,djql,sizeSubAp1,sizeSubAp2,cnh_hl,hl,l0h_hl,altlgs1,altlgs2,loworder=)
/* DOCUMENT
   
     <dipl> and <djql>   : arrays (nsubap,2) that contain the subaperture
                       coordinates projected onto altitude layer.
     <s1> and <s2>     : subap size of WFS in pupil plane
     <cnh_hl>              : scalar
     <l0h_hl>              : scalar, log10 of L0
     <altlgs1> and <altlgs2> : scalars, altitude of GUIDE STAR of wfs
     <hl>        : scalar, altitude of LAYER
     
   SEE ALSO:
*/
{
  // dipl and djql are the subaperture coordinates projected onto altitude layer
  dipl_x = dipl(,1);
  dipl_y = dipl(,2);
  djql_x = djql(,1);
  djql_y = djql(,2);
  
  // matrix of distances in x and y for all subapertures couples between 2 WFSs.
  Dijpql_x = (dipl_x(,-)-djql_x(-,));
  Dijpql_y = (dipl_y(,-)-djql_y(-,));

  sizeSubAp1_pup = sizeSubAp1;
  sizeSubAp2_pup = sizeSubAp2;
  
  // ..... CONE EFFECT IMPLEMENTATION ..... //
  if( altlgs1 ){
    sizeSubAp1 = sizeSubAp1*(altlgs1-hl)/altlgs1;
  }
  if( altlgs2 ){
    sizeSubAp2 = sizeSubAp2*(altlgs2-hl)/altlgs2;
  }

  nValidSubAp1 = dimsof(Dijpql_x)(2); //number of sub-apertures of 1st WFS
  nValidSubAp2 = dimsof(Dijpql_x)(3); //number of sub-apertures of 2nd WFS

  // test if the altitude layers is higher than the LGS altitude 
  if(sizeSubAp1<=0 || sizeSubAp2 <=0)
    return array(0., 2*nValidSubAp1, 2*nValidSubAp2);

  d1 = sizeSubAp1/2.;
  d2 = sizeSubAp2/2.;
  ac = d1-d2;
  ad = d1+d2;
  bc = -ad;
  bd = -ac;
  
  if(loworder){
    x0 = dm.pitch;
    if(sizeSubAp1 != sizeSubAp2){
      // Caa x-x
      caa_xx = -dphi_lowpass(Dijpql_x+ac,Dijpql_y,l0h_hl,x0) + dphi_lowpass(Dijpql_x+ad,Dijpql_y,l0h_hl,x0) + dphi_lowpass(Dijpql_x+bc,Dijpql_y,l0h_hl,x0) - dphi_lowpass(Dijpql_x+bd,Dijpql_y,l0h_hl,x0);
      // Caa y-y
      caa_yy = -dphi_lowpass(Dijpql_x,Dijpql_y+ac,l0h_hl,x0) + dphi_lowpass(Dijpql_x,Dijpql_y+ad,l0h_hl,x0) + dphi_lowpass(Dijpql_x,Dijpql_y+bc,l0h_hl,x0) - dphi_lowpass(Dijpql_x,Dijpql_y+bd,l0h_hl,x0);
    }else{
      //in this case, ac = bd = 0
      // Caa x-x
      caa_xx = -2*dphi_lowpass(Dijpql_x,Dijpql_y,l0h_hl,x0) + dphi_lowpass(Dijpql_x+ad,Dijpql_y,l0h_hl,x0) + dphi_lowpass(Dijpql_x+bc,Dijpql_y,l0h_hl,x0);
      // Caa y-y
      caa_yy = -2*dphi_lowpass(Dijpql_x,Dijpql_y,l0h_hl,x0) + dphi_lowpass(Dijpql_x,Dijpql_y+ad,l0h_hl,x0) + dphi_lowpass(Dijpql_x,Dijpql_y+bc,l0h_hl,x0);
    }
    // Caa x-y
    caa_xy = -dphi_lowpass(Dijpql_x+d1,Dijpql_y-d2,l0h_hl,x0) + dphi_lowpass(Dijpql_x+d1,Dijpql_y+d2,l0h_hl,x0) + dphi_lowpass(Dijpql_x-d1,Dijpql_y-d2,l0h_hl,x0) - dphi_lowpass(Dijpql_x-d1,Dijpql_y+d2,l0h_hl,x0);
  }else{
    if(sizeSubAp1 != sizeSubAp2){
      // Caa x-x
      caa_xx = -DPHI(Dijpql_x+ac,Dijpql_y,l0h_hl) + DPHI(Dijpql_x+ad,Dijpql_y,l0h_hl) + DPHI(Dijpql_x+bc,Dijpql_y,l0h_hl) - DPHI(Dijpql_x+bd,Dijpql_y,l0h_hl);
      // Caa y-y
      caa_yy = -DPHI(Dijpql_x,Dijpql_y+ac,l0h_hl) + DPHI(Dijpql_x,Dijpql_y+ad,l0h_hl) + DPHI(Dijpql_x,Dijpql_y+bc,l0h_hl) - DPHI(Dijpql_x,Dijpql_y+bd,l0h_hl);
    }else{
      // Caa x-x
      caa_xx = -2*DPHI(Dijpql_x,Dijpql_y,l0h_hl) + DPHI(Dijpql_x+ad,Dijpql_y,l0h_hl) + DPHI(Dijpql_x+bc,Dijpql_y,l0h_hl);
      // Caa y-y
      caa_yy = -2*DPHI(Dijpql_x,Dijpql_y+ac,l0h_hl) + DPHI(Dijpql_x,Dijpql_y+ad,l0h_hl) + DPHI(Dijpql_x,Dijpql_y+bc,l0h_hl);
    }
    // Caa x-y
    caa_xy = -DPHI(Dijpql_x+d1,Dijpql_y-d2,l0h_hl) + DPHI(Dijpql_x+d1,Dijpql_y+d2,l0h_hl) + DPHI(Dijpql_x-d1,Dijpql_y-d2,l0h_hl) - DPHI(Dijpql_x-d1,Dijpql_y+d2,l0h_hl);
  }
  
  caa = array(0.,2*nValidSubAp1,2*nValidSubAp2);
  caa(1:nValidSubAp1,1:nValidSubAp2) = caa_xx;
  caa(1:nValidSubAp1, 1+nValidSubAp2:nValidSubAp2*2) = caa_xy;
  caa(1+nValidSubAp1:nValidSubAp1*2,1:nValidSubAp2) = caa_xy;
  caa(1+nValidSubAp1:nValidSubAp1*2,1+nValidSubAp2:nValidSubAp2*2) = caa_yy;

  //units
  caa *=  0.5*(radian2arcsec*atm.lambda/2/pi)^2*(abs(cnh_hl))/sizeSubAp1_pup/sizeSubAp2_pup;
 
  return caa;
}

func covMat1layer(nbWfs,cnh_hl,hl,l0h_hl,xshift,yshift,magnification,theta,loworder=)
/* DOCUMENT calcCaa_3(nbWfs, L0, r0, altitude)
     
   <nbWfs>      : scalar, number of wfss
   <L0>         : scalar, log10(outer scale in meters)
   <r0>         : array 1D, lenght = # of layers
   <altitude>       : array 1D, lenght = # of layers
   
   
   SEE ALSO:
 */
{
  
  // reservation memoire
  covLearn_hl = array(0.,rtc.nSlopes,rtc.nSlopes);
 
  //managing the outerscale
  l0h_hl = abs(l0h_hl);

  if(l0h_hl >=100.){
    l0h_hl = 100.;
  }
  
  for(p=1; p<=nbWfs; p++) {       // loop 1 on WFSs
    // Compute coordinates of WFS #i on layer of h=<altitude>
    dipl = computeDipl(wfs(p).nLenslet,tel.diam,hl,wfs(p).x/radian2arcsec,wfs(p).y/radian2arcsec,wfs(p).lgsH,xshift(p),yshift(p),magnification(p),theta(p),obs=tel.obs);
    
    for(q=p; q<=nbWfs; q++) {       // boucle 2 sur les asos
      djql = computeDipl(wfs(q).nLenslet,tel.diam,hl,wfs(q).x/radian2arcsec, wfs(q).y/radian2arcsec,wfs(q).lgsH,xshift(q),yshift(q),magnification(q),theta(q),obs=tel.obs);


      if(wfs(p).type!=0 && wfs(q).type!=0 ) {   // if both WFS are valid ones
       
        // subap size in the pupil plane
        sspSize1 = tel.diam/wfs(p).nLenslet;
        sspSize2 = tel.diam/wfs(q).nLenslet;
        
        tmp = dij2covMat(dipl,djql,sspSize1,sspSize2,cnh_hl,hl,l0h_hl,wfs(p).lgsH,wfs(q).lgsH,loworder=loworder);
        
        rp = tslrange(p);
        rq = tslrange(q);
        covLearn_hl(rp, rq) = tmp;

      }
    } // end loop WFS 2
  }   // end loop WFS 1

  return covLearn_hl;
}

func covMatAllLayers(nbWfs,cnh,alt,l0h,xshift,yshift,magnification,theta,filtreTilt,loworder=,verb=)
/* DOCUMENT calcMatCov_explicit(nbWfs, L0, r0, alti, xaso, yaso, XPup, YPup, diamPup, thetaML)
   
   <nbWfs>      : scalar, number of wfss
   <L0>         : array 1D, scalar, log10(outer scale in meters), lenght = # of layers
   <r0>         : array 1D, lenght = # of layers
   <alti>       : array 1D, lenght = # of layers
   <filtreTilt> : tip-tilt is filtered on LGS when set to 1. Default=1.
   
   SEE ALSO:
   */
{
  local covLearn;
  
  covLearn = array(0.,rtc.nSlopes,rtc.nSlopes);

  // loop on turbulent layers
  nl = learn.nl;
  for(l=1; l<=nl; l++) {
    if(cnh(l) !=0){
      if(verb){
        write,format="\rGenerating the %d th layer...",l;
      }
      covLearn += covMat1layer(nbWfs,cnh(l),alt(l),l0h(l),xshift,yshift,magnification,theta,loworder=loworder);
    }
  }

  // Filtering of tiptilt according to the WFS-type : 0, 2 or 3 (no wfs, lgs, TTonly)

  if( is_void(filtreTilt ) ) filtreTilt = 1;
  if( filtreTilt==1 ) {
    covLearn = handle_tilt_from_wfstype(covLearn,ttr=learn.ttr);
  }else{
    if(!singleWFS)
      covLearn =  duplicateCovar(covLearn);
  }
  return covLearn;
  
}

func covMatModel(learn, fitEstim,loworder=,verb=)
/* DOCUMENT 
     
     
   SEE ALSO:
 */
{
  local cnh,alt,l0h,tracking,xshift,yshift,magnification,theta;
  
  // get all the variables that were encapsulated inside fitEstim
  unpackcoeffs,fitEstim, cnh, alt, l0h, tracking, xshift, yshift,magnification,theta;
  
  covLearn = covMatAllLayers(rtc.nWfs,cnh,alt,l0h,xshift,yshift,magnification,theta,filtreTilt,loworder=loworder);

   // ................... Tracking perturbation ........................
  if(!learn.ttr){
     covLearn  = trackingMatCov(tracking,  covLearn );
  }
  
  // 2-step case
  // this matrix is useful only when doing a "2-step Learn", for removing the ground layer
  if(learn.runLearnTwoSteps){
    P = *learn.transformationMatrix;
     covLearn  = (P(,+) *  covLearn (,+))(,+) *P(+,);
  }
 
  // ................... Noise options (diagonal or not)  .......................
  if(learn.diagonal==1) {
    takesDiag( covLearn ) = 0.;
  }

  return  covLearn ;
}

func trackingMatCov(tracking,covLearn)
/* DOCUMENT trackingMatCov(tracking, fullMat)
     Creates a checkerboard on top of the covariance matrix, with
     <x^2>=tracking(1) on all the x/x covariance areas,
     <y^2>=tracking(2) on all the y/y covariance areas,
     <xy> =tracking(3) on all the x/y and y/x covariance areas.
   SEE ALSO:
 */
{
  fullMat = covLearn;
  for(p=1;p<=rtc.nWfs;p++) {
    for(q=p;q<=rtc.nWfs;q++) {
      if( wfs(p).type!=2 && wfs(q).type!=2 ) {
        // ranges for X-slopes
        xrp = tslindex(p):(tslindex(p)+wfs(p).nValidSubAp-1);
        xrq = tslindex(q):(tslindex(q)+wfs(q).nValidSubAp-1);
        // ranges for Y-slopes
        yrp = tslindex(p)+wfs(p).nValidSubAp:tslindex(-p);
        yrq = tslindex(q)+wfs(q).nValidSubAp:tslindex(-q);
        // addition of constants
        fullMat( xrp, xrq ) += tracking(1);
        fullMat( yrp, yrq ) += tracking(2);
        fullMat( xrp, yrq ) += tracking(3);
        fullMat( yrp, xrq ) += tracking(3);
        
        if(q!=p){
          fullMat( xrq, xrp ) += tracking(1);
          fullMat( yrq, yrp ) += tracking(2);
          fullMat( xrq, yrp ) += tracking(3);
          fullMat( yrq, xrp ) += tracking(3);
        }
      }
    }
  }
 
  
  return fullMat;
}

/*
____  _                      ____  _                   _                  
|  _ \| |__   __ _ ___  ___  / ___|| |_ _ __ _   _  ___| |_ _   _ _ __ ___ 
| |_) | '_ \ / _` / __|/ _ \ \___ \| __| '__| | | |/ __| __| | | | '__/ _ \
|  __/| | | | (_| \__ \  __/  ___) | |_| |  | |_| | (__| |_| |_| | | |  __/
|_|   |_| |_|\__,_|___/\___| |____/ \__|_|   \__,_|\___|\__|\__,_|_|  \___|
                                                                           
 _____                 _   _                 
|  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
| |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
|  _|| |_| | | | | (__| |_| | (_) | | | \__ \
|_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
                                             
*/

func DPHI(x,y,L0,model=)
{
  if(is_void(model)){model = "Von-Karman";}
  if(L0 >= 1e3){model = "Kolmogorov";}
  
  r = abs(x,y);
  
  if(model == "Von-Karman"){
    return rodconan(r, L0);
  }else if(model == "Kolmogorov"){
    return 6.88 * r^(5/3.);
  }
}



func rodconan(r,L0,k=)
/* DOCUMENT rodconan(r,L0,k=)
     The phase structure function is computed from the expression
     Dphi(r) = k1  * L0^(5./3) * (k2 - (2.pi.r/L0)^5/6 K_{5/6}(2.pi.r/L0))

     For small r, the expression is computed from a development of
     K_5/6 near 0. The value of k2 is not used, as this same value
     appears in the series and cancels with k2.
     For large r, the expression is taken from an asymptotic form.

     
     Computes the phase structure function for a separation <r>.  The
     r0 is not taken into account : the final result of rodconan(r,L0)
     has to be scaled with r0^-5/3, with r0 expressed in same units as
     <r>, to get the right value.

     dphi_r = rodconan(r, L0) * r0^(5./3);
   
   SEE ALSO:
 */
{
  local k;
  // k1 is the value of :
  k1 = 2*gamma_R(11./6)*2^(-5./6)*pi^(-8./3)*(24*gamma_R(6./5)/5.)^(5./6);
  k2 = gamma_R(5./6)*2^(-1./6);
  
  dprf0 = (2*pi/L0)*r;

  res = r;
  Xlim = 0.75*2*pi;
  largeX = dprf0>Xlim;
  ilarge = where(largeX);
  ismall = where(!largeX);
  if( is_array(ilarge) ) {
    res(ilarge) = asymp_macdo(dprf0(ilarge));
  }
  if( is_array(ismall) ) {
    res(ismall) = -macdo_x56(dprf0(ismall), k=k);
  }
  return (k1 * L0^(5./3)) * res; 
}

func macdo_x56(x,k=)
/* DOCUMENT  macdo_x56(x)
   
     Computation of the function
     f(x) =*K_{5/6}(x)
     using a series for the esimation of K_{5/6}, taken from Rod Conan thesis :
     K_a(x)=1/2 \sum_{n=0}^\infty \frac{(-1)^n}{n!}
     \left(\Gamma(-n-a) (x/2)^{2n+a} + \Gamma(-n+a) (x/2)^{2n-a} \right) ,
     with a = 5/6.

     Setting x22 = (x/2)^2, setting uda = (1/2)^a, and multiplying by x^a,
     this becomes :
     x^a * Ka(x) = 0.5 $ -1^n / n! [ G(-n-a).uda x22^(n+a) + G(-n+a)/uda x22^n ]
     Then we use the following recurrence formulae on the following quantities :
     G(-(n+1)-a) = G(-n-a) / -a-n-1
     G(-(n+1)+a) = G(-n+a) /  a-n-1
     (n+1)! = n! * (n+1)
     x22^(n+1) = x22^n * x22
     and at each iteration on n, one will use the values already computed at step (n-1).
     The values of G(a) and G(-a) are hardcoded instead of being computed.

     The first term of the series has also been skipped, as it
     vanishes with another term in the expression of Dphi.
     
   SEE ALSO:
 */
{
  a = 5./6.;
  if( is_void(k) ) k=10;
  fn = 1.;                             // initialisation factorielle 0!=1
  x2a = x^(2.*a);
  x22 = x*x/4.;                        //  (x/2)^2
  x2n = 0.5;                           // init (1/2) * x^0
  Ga  = gamma_R(a) / (1/2.)^a;       // Gamma(a) / (1/2)^a
  Gma = gamma_R(-a) * (1/2.)^a;       // Gamma(-a) * (1/2.)^a
  s = array(0.0, dimsof(x));
  for(n=0; n<=k; n++) {
    dd = Gma * x2a;
    if( n )
      dd += Ga;
    dd *= x2n;
    dd /= fn;
    // addition to s, with multiplication by (-1)^n
    if( n%2 ) s -= dd;
    else      s += dd;
    // prepare recurrence iteration for next step
    if( n<k ) {
      fn *= n+1;     // factorial
      Gma /= -a-n-1; // gamma function
      Ga /= a-n-1;   // idem
      x2n *= x22;    // x^n
    }
  }
  return s;
}



func asymp_macdo(x)
/* DOCUMENT asymp_macdo(x)

     Computes a term involved in the computation of the phase struct
     function with a finite outer scale according to the Von-Karman
     model. The term involves the MacDonald function (modified bessel
     function of second kind) K_{5/6}(x), and the algorithm uses the
     asymptotic form for x ~ infinity.
     Warnings :
         - This function makes a floating point interrupt for x=0
     and should not be used in this case.
         - Works only for x>0.
     
   SEE ALSO:
 */
{
  k2 = gamma_R(5./6)*2^(-1./6)
  k3 = sqrt(pi/2);
  a1 = 2./9;
  a2 = -7/81.;  //  -7/81
  a3 = 175./2187;

  x_1 = 1./x;
  res = k2 - k3*exp(-x)*x^(1/3.)*(1.0 + x_1*(a1 + x_1*(a2 + x_1*a3)));
  return res;
}


func dphi_lowpass(x,y,L0,x0)
{
  r = abs(x,y);
  
  if(L0 >= 1e4){
    //model Kolmogorov
    return 6.1554764445358891578 * r^(5/3.) * Ij0t83(r*pi/x0);
  }else{
    //model Von-Karman
    return DPHI(x,y,L0) - dphi_highpass(r,x0);
  }
}
func dphi_highpass(r,x0)
{
  //write,format="%.20g\n",(2*(2*pi)^(8/3.)*0.0228956) -> 6.1554764445358891578
  return (r^(5./3.)) * (1.1183343328701949 - Ij0t83(r*pi/x0)) * 6.1554764445358891578;
}


func Ij0t83(x)
/* DOCUMENT 
     
   Calcul de l'integrale tabulee
    x
   $ t^(-8/3) (1-bessel_j(0,t)) dt
    0

    Pres de 0, le resultat s'approxime par (3/4.)*x^(1./3)*(1-x^2/112.+...)
    
   SEE ALSO:

*/
{
  if( numberof(dimsof(x))==1 ) {   // x is a scalar
    return Ij0t83([x])(1);
  }
  y = double(x);
  largeX = x>exp(-3.0);
  nn = where( largeX );
  if( is_array(nn) ) {
    y(nn) = prepare_large_Ij0t83(x(nn));
  }
  nn = where( !largeX );
  if( is_array(nn) ) {
    y(nn) = prepare_small_Ij0t83(x(nn));
  }
  return y;
}


func prepare_small_Ij0t83(x)
{
  return 0.75*x^(1./3)*(1-x^2/112.);
}

func prepare_large_Ij0t83(xp)
{
 
  n = 100;
  t = span(-4,10,n);   // u = exp(t);
  dt = (t(0)-t(1))/(n-1) ;
  A = prepare_small_Ij0t83( exp(-4.0) );   // integral from -inf to exp(-4.0)
  XLARGE = exp(t);
  y = exp(-t*(5/3.))*unMoinsJ0(XLARGE);  // u^(-8./3) * (1-J0(u)) du  avec du=exp(t).dt  
  YLARGE = y(zcen)(cum) * dt + A;

  return interp(YLARGE,XLARGE,xp);
}

func unMoinsJ0(x)
{
  if( numberof(dimsof(x))==1 ) {   // x is a scalar
    if( x<0.1 ) {
      // J0(x) = x^2/4 - x^4/64 + ...
      //       = (x/2)^2  ( 1 - (x/2)^2 / 4 + ... ) to minimize computation errors
      x22 = (x/2.)^2; 
      return (1-x22/4.)*x22;
    } else {
      // classique
      return (1-bessj0(x));
    }
  } else {  // x is a vector
    y = double(x);
    for(i=1; i<=numberof(x); i++)
      y(i) = unMoinsJ0(x(i));
    return y;
  }
}

func gamma_R(x)
/* DOCUMENT gamma(x) : fonction gamma generalisee a R
   utilisation de la relation :
   gamma(x)*gamma(-x)=-pi/(x*sin(pi*x))
*/
{
  require,"gamma.i";
  
  //on impose gamma(1)=1;
  gamma_1=1.0;
  
  if (is_scalar(x))
    {
      a=[x];
    }
  else
    {
      a=x;
    }

  a*=1.0;
  result=a;
  
  //on fait le test sur le domaine de a
  idx=where((a>0.0)&(a<1.0));
  if (numberof(idx)!=0)
    {
      z=a(idx);
      result(idx)=exp(ln_gamma(z));
    }
  
  //on fait le test pour a==1
  idx=where((a==1.0));
  if (numberof(idx)!=0)
    result(idx)=gamma_1;
  
  //on fait le test pour a>=1
  idx=where(a>=1.0);
  if (numberof(idx)!=0)
    {
      z=a(idx);
      result(idx)=exp(lngamma(z));
    }
  
  //on fait le test pour -1<z<0
  idx=where((a>-1.0)&(a<0.0));
  if (numberof(idx)!=0)
    {
      z=-1.0*a(idx);
      temp=-1.0*pi;
      temp/=z;
      temp/=sin(pi*z);
      temp/=exp(ln_gamma(z));
      result(idx)=temp;
    }
  
  //on fait le test pour z<-1
  idx=where(a<-1.0);
  if (numberof(idx)!=0)
    {
      z=-1.0*a(idx);
      temp=-1.0*pi;
      temp/=z;
      temp/=sin(pi*z);
      temp/=exp(lngamma(z));
      result(idx)=temp;
    }
  
  if (is_scalar(x))
    return(result(1));
  else
    return(result);
}

/*
 __  __       _        _               
|  \/  | __ _| |_ _ __(_) ___ ___  ___ 
| |\/| |/ _` | __| '__| |/ __/ _ \/ __|
| |  | | (_| | |_| |  | | (_|  __/\__ \
|_|  |_|\__,_|\__|_|  |_|\___\___||___/
                                       
 __  __                                                   _   
|  \/  | __ _ _ __   __ _  __ _  ___ _ __ ___   ___ _ __ | |_ 
| |\/| |/ _` | '_ \ / _` |/ _` |/ _ \ '_ ` _ \ / _ \ '_ \| __|
| |  | | (_| | | | | (_| | (_| |  __/ | | | | |  __/ | | | |_ 
|_|  |_|\__,_|_| |_|\__,_|\__, |\___|_| |_| |_|\___|_| |_|\__|
                          |___/                               
*/


func Filt_Zer(modefilt,filt_ngs=,filt_lgs=,filt_ttgs=)
/* DOCUMENT MatFiltZer = Filt_Zer(modefilt,filt_ngs=,filt_lgs=,filt_ttgs=);
   Cette fonction genere une matrice qui filtre les modefilt premiers modes de Zernike des donnees.
   Ex: si on veut filtrer le tip-tilt, on calcule MatFiltZer = Filt_Zer(2) et les donnees sans TT se deduise par l'operation st = MatFiltZer(,+)*s(+,).
   Si filt_xgs = 1, le filtrage se fait sur les aso des type x (N, L ou TT)
   Si filt_xgs = 0, pas de filtrage (identite sur la partie de la matrice a considerer)
   Si filt_xgs = -1, mise a 0 dans la matrice


   Olivier Martin
 */

{
  
  MatFiltZer = identite(rtc.nSlopes);
  if(!filt_ngs) filt_ngs=0;
  if(!filt_lgs) filt_lgs=0;
  if(!filt_ttgs) filt_ttgs=0;
  
  if(modefilt!=0){
    
    for(i=1; i<=rtc.nWfs; i++) {
      
      if((wfs(i).type==1 && filt_ngs==1) || (wfs(i).type==2 && filt_lgs==1)  || (wfs(i).type==3 && filt_ttgs==1) ) {
        mrz = *sys.slopesToZernikeMatrix;
        mrz(modefilt+1:,) = 0;

        zmi = (*sys.zernikeToSlopesMatrix)(,1:modefilt);
        mrz = invgen(zmi,0, t=1);
        nslopes_i = wfs(i).nValidMeas;
        MatFiltZer(slrange(i),slrange(i)) = identite(nslopes_i) - (zmi)(,+) * mrz(+,);
      }
      // si ASO desactive on vire tout
      if(wfs(i).type==0 || (wfs(i).type==1 && filt_ngs==-1) || (wfs(i).type==2 && filt_lgs==-1)  || (wfs(i).type==3 && filt_ttgs==-1))
        MatFiltZer(slrange(i),slrange(i)) = 0;
    }
  }

  return MatFiltZer;
}


func Keep_Zer(modefilt,calc_ngs=,calc_lgs=,calc_ttgs=)
/* DOCUMENT MatKeepZer = Keep_Zer(modefilt);
   
   Cette fonction genere une matrice qui calcule les modefilt premiers modes de Zernike des donnees.
   Ex: si on veut garder seulement le tip-tilt, on calcule MatKeepZer = Keep_Zer(2) et le TT se deduit par l'operation st = MatKeepZer(,+)*s(+,).
   Si calc_ngs=1, le calcul se fait sur les ngs
   Si calc_lgs=1, le calcul se fait sur les lgs
   Si les deux valent 1, le calcul se fait sur tous les asos NGS ou LGS

   Olivier Martin

 */

{

  if(modefilt!=0){
    
    if(!calc_ngs) calc_ngs=0;
    if(!calc_lgs) calc_lgs=0;
    if(!calc_ttgs) calc_ttgs=0;

    MatKeepZer = identite(rtc.nSlopes) - Filt_Zer(modefilt,filt_ngs=calc_ngs,filt_lgs=calc_lgs,filt_ttgs=calc_ttgs);
    return MatKeepZer;
    
  }
  
  else return "modefilt is equal to 0";  

}

func Calc_CommonMat(avg_ngs=,avg_lgs=,avg_ttgs=)
/* DOCUMENT MatCompGL = Calc_CommonMat()
   Cette fonction genere une matrice qui prend la moyenne des pentes sur les Asos. Si sgl = MatCompGL(,+)*s(+,),
   sgl represente la couche au sol:
   sgl(i,) = sum_{k=1}^{k=rtc.nbWfs} s(i+k*rtc.Nslopes,)/rtc.nbWfs

   Si avg_xgs=1, on moyenne sur les WFS de type x (N,L ou TT)  

   Olivier Martin
 */

{

  // matrice qui permet de prendre la moyenne sur les ASOs
  MatCompGL = array(0.,rtc.nSlopes,rtc.nSlopes);
  if(!avg_ngs) avg_ngs=0;
  if(!avg_lgs) avg_lgs=0;
  if(!avg_ttgs) avg_ttgs=0;
  
  activeWfs = (wfs.type==1 & avg_ngs==1) | (wfs.type==2 & avg_lgs==1) | (wfs.type==3 & avg_ttgs==1);
  
  nbActiveWfs = 0;
  for(i=1;i<=rtc.nWfs;i++) {
    if( activeWfs(i) ) {
      nbActiveWfs++;
      for(j=1;j<=rtc.nWfs;j++) {
        if( activeWfs(j) ) {
          MatCompGL(slrange(i),slrange(j)) = MatCompGL(slrange(j),slrange(i)) = identite(wfs(i).nValidMeas);
        }
      }
    }
  }

  if(nbActiveWfs==0)
    return MatCompGL;

  MatCompGL /= nbActiveWfs;  
  return MatCompGL;
}




func ComputeTransMatrix2steps(modefilt,ttr=)
/* DOCUMENT P = ComputeTransMatrix2steps(modefilt)
   
  Cette fonction permet de virer la couche au sol en ne considerant
  pas les modefilt premiers modes de Zernike des LGS. Cette fonction
  est utilisee dans la procedure 2 steps et la matrice de covariance
  des couches en altitude se calcule par:
  caa_alt = (P(,+)*caa_tot(,+))(,+)*P(+,);

   SEE ALSO Calc_Zer_NGS, Calc_CommonMat, Filt_Zer

   Olivier Martin
 */


{
  if(!ttr){
    //1. On filtre le TT sur toutes les mesures
    MatKeepHO = Filt_Zer(modefilt,filt_ngs=1,filt_lgs=1,filt_ttgs=-1);
    //2. On prend la moyenne des WFS sans TT
    MatCompGL_HO = Calc_CommonMat(avg_ngs=1,avg_lgs=1,avg_ttgs=0);
    //3. On garde la partie TT des WFS de type 1 et 3
    MatKeepTT = Keep_Zer(modefilt,calc_ngs=1,calc_lgs=0,calc_ttgs=1);
    //4.On prend la moyenne sur les hauts ordres
    MatCompGL_TT = Calc_CommonMat(avg_ngs=1,avg_lgs=0,avg_ttgs = 1);  
    //4. On vire la couche au sol des mesures
    P =  identite(rtc.nSlopes) - (MatCompGL_TT(,+)*MatKeepTT(,+) + MatCompGL_HO(,+)*MatKeepHO(,+));
  }else{
    //We only average slopes that are all tip-tilt removed
    MatKeepHO = Filt_Zer(modefilt,filt_ngs=1,filt_lgs=1,filt_ttgs=-1);
    MatCompGL_HO = Calc_CommonMat(avg_ngs=1,avg_lgs=1,avg_ttgs=0);
    P =  identite(rtc.nSlopes) -  MatCompGL_HO(,+)*MatKeepHO(,+);
  }
  
  return P;

}





func duplicateCovar( caa_input )
{

  for(i=1; i<=rtc.nWfs; i++) {
    for(j=i; j<=rtc.nWfs; j++) {
      
      ri = tslrange(i);
      rj = tslrange(j);
      
      if( wfs(i).type==0 || wfs(j).type==0 ) {
        caa_input(ri,rj) = 0.0;
        caa_input(rj,ri) = 0.0;
      } else {
        // Duplicating the matrix on the other side of the diagonal
        if( i!=j )
          caa_input(rj,ri) = transpose(caa_input(ri,rj));
      }
    }
    
  }
  return caa_input;
}


func handle_tilt_from_wfstype( caa_input, ttr= )
/* DOCUMENT handle_tilt_from_wfstype( Caa )

   Tilt filtering of a covariance matrix of ALL the wfs (incl. TS).
   The tilt-filtering process is adapted to the WFS type, either LGS
   (tilt removed) or tilt-only (HO removed), or devalidated wfs (all
   set to 0).

   WARNING : the matrix caa_full MUST be symmetric !!
   
   SEE ALSO:
 */
{
  // this line is required if we don't want to modify the input array
  caa_full = caa_input;
  if(is_void(ttr)) ttr = 0;

  KeepTT = filt_TTvec(wfs(rtc.its).nValidMeas);
  FiltTT = calc_TTMatFilt_1(wfs(rtc.its).nValidMeas);
  
  for(i=1; i<=rtc.nWfs; i++) {
    for(j=i; j<=rtc.nWfs; j++) {
      
      ri = tslrange(i);
      rj = tslrange(j);
      
      tmp = caa_full(ri,rj);
      if(!ttr){
        if(wfs(i).type==1) {
          tmpi = tmp;
        } else if(wfs(i).type==2) {
          tmpi = FiltTT(,+)*tmp(+,);
        } else if(wfs(i).type==3) {
          tmpi =  KeepTT(,+)*tmp(+,);
        } else {
          tmpi = tmp * 0.0;
        }
        if(wfs(j).type==1) {
          caa_full(ri,rj) = tmpi;
        } else if(wfs(j).type==2) {
          caa_full(ri,rj) = tmpi(,+)*FiltTT(+,);
        } else if(wfs(j).type==3) {
          caa_full(ri,rj) = tmpi(,+) *KeepTT(+,);
        } else {
          caa_full(ri,rj) = tmpi * 0.0;
        }
      }else{
        tmpi = FiltTT(,+)*tmp(+,);
        caa_full(ri,rj) = tmpi(,+) * FiltTT(+,);
      }
      
      // Duplicating the matrix on the other side of the diagonal
      if( i!=j )
        caa_full(rj,ri) = transpose(caa_full(ri,rj));
    }
    
  }
  return caa_full;
}

func calc_TTMatFilt_1(nbslopes)
/* DOCUMENT returns a matrix <p1> that filters the average TT on the vector of slopes <vec>
   with nbslopes = x+y slopes
   
   The filtered vector can be computed using:
   vec_filt = p1(,+)*vec(+)
     
   SEE ALSO: filt_TTvec
 */
{
  p = array( -2./nbslopes, nbslopes, nbslopes);
  p(1:nbslopes/2, 1+nbslopes/2:) = 0.0;
  p(1+nbslopes/2:,1:nbslopes/2) = 0.0;
  p(*)(1::nbslopes+1) += 1;
  return p;
}

func filt_TTvec(nbslopes)
/* DOCUMENT returns a matrix p that averages the slopes over x and y. 

   The filtered vector can be computed using:
   vec_filt = p(,+)*vec(+)
   
   SEE ALSO: calc_TTMatFilt_1
 */
{
  p = array( 2./nbslopes, nbslopes, nbslopes);
  p(1:nbslopes/2, 1+nbslopes/2:) = 0.0;
  p(1+nbslopes/2:,1:nbslopes/2) = 0.0;
  return p;
}

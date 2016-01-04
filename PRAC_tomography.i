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
  Cnn = array(0.,sum(data.wfs.nssp),sum(data.wfs.nssp));
  takesDiag(Cnn) = *data.turbu.varNoise;
  
  if(obsmode == "MOAO"){
    Cnn =  handle_tilt_from_wfstype(Cnn);
    data.covmatrix.covNoise = &Cnn;
    if(verb){
      write,"\rComputing covariance matrices of measurements...";
    }
    
    //computes the synthetic covariances matrices of WFS measurements
    Caa = computesCaaMatrix(data.nwfs,data.learn.cnh,data.learn.altitude,data.learn.l0h,data.wfs.x,data.wfs.y,verb=verb);

    Cpara2 = *data.covmatrix.covPara;
    Cpara2(norange(data.its),norange(data.its)) += Caa(norange(data.its),norange(data.its));
    
    if(verb){
      write,"\rComputing tomographic residual covariance matrix...";
    }
    
    Cee = computeCovSlopesError(Cpara2,*data.rtc.R);
    //managing the aliasing
    Caa(slrange(data.its),) = 0;
    Caa(,slrange(data.its)) = 0;
    Cee += computeCovSlopesError(Caa,*data.rtc.R);
    
    //manages the time-filtered off-axis noise contribution
    noisefact = filteringNoiseFactor(data.rtc.gain,data.rtc.delay*data.rtc.Fe,data.rtc.obsmode);
    data.budget.noise = computesTomoError(Cnn,*data.rtc.R);
    Cee += noisefact  * computeCovSlopesError(Cnn,*data.rtc.R);

  }else if(obsmode == "SCAO"){
    
    if(verb){
      write,"\rComputing the aliasing and noise injected into the loop..."
    }
    //noise propgated through the loop
    fgt = filteringNoiseFactor(data.rtc.gain,data.rtc.delay*data.rtc.Fe,data.rtc.obsmode);
    Cee =  fgt*Cnn(slrange(data.its),slrange(data.its));
    
    //computes the aliasing covariance matrix
    Caa = computesCaaMatrix(1,data.learn.cnh,data.learn.altitude,data.learn.l0h,data.wfs.x,data.wfs.y,singleWFS=data.its);
    
    //manages the time-filtering of the aliasing
    Nframes = data.nframes;
    nu = data.rtc.Fe*(indgen(Nframes)-1)/Nframes;
    hbf = abs(hbfScao(nu,data.rtc.Fe,data.rtc.delay,data.rtc.gain,data.rtc.BP))^2;
    aliasfact = sum(hbf)/numberof(hbf);
    Cee += Caa;

    if(verb){
      write,format="\rNoise filtering factor: %g\n",noisefact;
      write,format="\rAliasing filtering factor: %g\n",aliasfact;
    }

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
  
  rr = 1:tslindex(1-data.nwfs);  // WFS 1 to last before TS
  tsr = tslrange(data.its);
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
    if(is_void(TScam)) TScam = data.its;
    YY = MatCov(slrange(TScam),slrange(TScam));
    XY = MatCov(norange(TScam),slrange(TScam));
    XX = MatCov(norange(TScam),norange(TScam));
          
    //Calcul erreur minimale avec la covariance du TS
    Cee = YY - XY(+,)*R(,+) - R(,+)*XY(+,) + R(,+)*(XX(,+)*R(,+))(+,);
  }else{
    Cee = MatCov;
  }
  
  //Derivation on the Zernike
  mrz = *data.wfs(TScam).mrz;
  Czz = mrz(,+)*(Cee(,+)*mrz(,+))(+,);
  nzz = dimsof(Czz)(0);
  // units ...
  units = ((1/206265.)*(0.5*(data.tel.diam)*1e9));
  diagonale = takesDiag(Czz)*units^2;
  
  if(full){
    return SQRT(diagonale);
  }else{
    return SQRT((sum(diagonale)));
  }
}

func computesCaaMatrix(nwfs,cnh,alt,l0h,posX,posY,singleWFS=,verb=)
{
  if(singleWFS){
    filtreTilt = 0;
  }
  nwfstmp = data.nwfs;
  data.nwfs = nwfs;
  diaginit = data.learn.diagonal;
  data.learn.diagonal = 0;
  data.learn.ttr=0;
  //computes the covariance matrix of the whole turbulence
  fitEstim = packcoeffs(data);
  
  Cturbu = covMatModel(data, fitEstim,loworder=0,singleWFS=singleWFS,verb=verb);
  data.covmatrix.covLearn = &Cturbu;

  //computes the covariance matrix of the low order modes only
  Cpara = covMatModel(data, fitEstim,loworder=1,singleWFS=singleWFS,verb=verb);
  data.covmatrix.covPara = &Cpara;
  
  data.nwfs = nwfstmp;
  data.learn.diagonal = diaginit;

  data.covmatrix.covAlias = &(Cturbu - Cpara);

  //saving tracking covariance matrix.
  covTracking = 0*Cturbu;
  if(data.rtc.obsmode == "MOAO"){
    covTracking = trackingMatCov(data.learn.tracking, covTracking);
  }else if(data.rtc.obsmode == "SCAO"){
    ns = dimsof(covTracking)(0);
    covTracking += data.learn.tracking(3);
    covTracking(1:ns/2,1:ns/2) = data.learn.tracking(1);
    covTracking(ns/2+1:,ns/2+1:) = data.learn.tracking(2);
  }
  data.covmatrix.covTracking = &covTracking;
    
  return *data.covmatrix.covAlias;
}

/*
  ____                           _              
 / ___| ___  ___  _ __ ___   ___| |_ _ __ _   _ 
| |  _ / _ \/ _ \| '_ ` _ \ / _ \ __| '__| | | |
| |_| |  __/ (_) | | | | | |  __/ |_| |  | |_| |
 \____|\___|\___/|_| |_| |_|\___|\__|_|   \__, |
                                          |___/ 
*/

func getValidSubapArray( nssp, rext, rint )
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
  
  x=span(-1,1,nssp+1)(zcen)(,-:1:nssp);
  y=transpose(x);
  r=sqrt(x*x+y*y);
  valid2d = (r<rext & r>=rint);
  valid = valid2d(*);
  return valid;
}

func compute_ri(nssp, diamtel, obs,xshift,yshift)
/* DOCUMENT vectmesFab(nssp, diamtel, obs)
   
     
   SEE ALSO:
 */
{
  rtel = diamtel/2.;
  x = span(-rtel, rtel, nssp+1)(zcen)(,-:1:nssp);
  y = transpose(x);
  nn = where(getValidSubapArray( nssp, 1.0, obs ));
  // taking pupil offsets (misregistration deviations) into account
  x += xshift;
  y += yshift;
  return [x(*)(nn),y(*)(nn)];
}


func compute_dipl(np, D, h, x, y, Hlgs, xshift,yshift,obs=)
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
  
  dipl = compute_ri(np,D, obs, xshift,yshift)(..,-:1:n);
  
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

func Dij2covMat(dipl, djql, s1, s2, cnh_hl, hl, l0h_hl, altlgs1, altlgs2,p,q,loworder=,deriv=)
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

  s1_pup = s1;
  s2_pup = s2;
  // rescale the subaperture size the projected size of a subaperture at a given altitude layer
  if( altlgs1 ){
    s1 = s1*(altlgs1-hl)/altlgs1;
  }
  if( altlgs2 ){
    s2 = s2*(altlgs2-hl)/altlgs2;
  }

  nssp1 = dimsof(Dijpql_x)(2); //number of sub-apertures of 1st WFS
  nssp2 = dimsof(Dijpql_x)(3); //number of sub-apertures of 2nd WFS
  // test if the altitude layers is higher than the LGS altitude 
  if(s1<=0 || s2 <=0) return array(0., 2*nssp1, 2*nssp2);

  d1 = s1/2;
  d2 = s2/2;
  ac = d1-d2;
  ad = d1+d2;
  bc = -ad;
  bd = -ac;
  
  
  if(loworder){
    x0 = data.tel.diam/data.wfs(data.its).sX;
    if(s1 != s2){
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
    if(!deriv){
      if(s1 != s2){
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

    }else{
      if(s1 != s2){
        // Caa x-x
        caa_xx = -gradDPHI(Dijpql_x+ac,Dijpql_y,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x+ad,Dijpql_y,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x+bc,Dijpql_y,cnh_hl,l0h_hl,p,q,deriv=deriv) - gradDPHI(Dijpql_x+bd,Dijpql_y,cnh_hl,l0h_hl,p,q,deriv=deriv);
        // Caa y-y
        caa_yy = -gradDPHI(Dijpql_x,Dijpql_y+ac,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x,Dijpql_y+ad,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x,Dijpql_y+bc,cnh_hl,l0h_hl,p,q,deriv=deriv) - gradDPHI(Dijpql_x,Dijpql_y+bd,cnh_hl,l0h_hl,p,q,deriv=deriv);
      }else{
        // Caa x-x
        caa_xx = -2*gradDPHI(Dijpql_x,Dijpql_y,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x+ad,Dijpql_y,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x+bc,Dijpql_y,cnh_hl,l0h_hl,p,q,deriv=deriv);
        // Caa y-y
        caa_yy = -2*gradDPHI(Dijpql_x,Dijpql_y+ac,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x,Dijpql_y+ad,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x,Dijpql_y+bc,cnh_hl,l0h_hl,p,q,deriv=deriv);
      }
      // Caa x-y
      caa_xy = -gradDPHI(Dijpql_x+d1,Dijpql_y-d2,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x+d1,Dijpql_y+d2,cnh_hl,l0h_hl,p,q,deriv=deriv) + gradDPHI(Dijpql_x-d1,Dijpql_y-d2,cnh_hl,l0h_hl,p,q,deriv=deriv) - gradDPHI(Dijpql_x-d1,Dijpql_y+d2,cnh_hl,l0h_hl,p,q,deriv=deriv); 
    }
  }
  
  caa = array(0., 2*nssp1, 2*nssp2);
  caa(1:nssp1,1:nssp2) = caa_xx;
  caa(1:nssp1, 1+nssp2:nssp2*2) = caa_xy;
  caa(1+nssp1:nssp1*2,1:nssp2) = caa_xy;
  caa(1+nssp1:nssp1*2,1+nssp2:nssp2*2) = caa_yy;

  //units
  caa *=  0.5*(206265*0.5e-6/2/pi)*(206265*0.5e-6/2/pi)*(abs(cnh_hl))/s1_pup/s2_pup;
 
  return caa;
}

func covMat1layer(nbWfs, l,cnh_hl, hl, l0h_hl,xshift,yshift,loworder=,singleWFS=,deriv=)
/* DOCUMENT calcCaa_3(nbWfs, L0, r0, altitude)
     
   <nbWfs>      : scalar, number of wfss
   <L0>         : scalar, log10(outer scale in meters)
   <r0>         : array 1D, lenght = # of layers
   <altitude>       : array 1D, lenght = # of layers
   
   
   SEE ALSO:
 */
{
  
  // reservation memoire
  if(!singleWFS){
    Cmod_hl = array(0.,sum(data.wfs.nssp), sum(data.wfs.nssp));
  }else{
    Cmod_hl = array(0.,data.wfs(singleWFS).nssp, sum(data.wfs(singleWFS).nssp));
  }

  //managing the outerscale
  l0h_hl = abs(l0h_hl);

  if(l0h_hl >=100.){
    l0h_hl = 100.;
  }
  
  for(p=1; p<=nbWfs; p++) {       // loop 1 on WFSs
    // Compute coordinates of WFS #i on layer of h=<altitude>
    if(singleWFS){m = singleWFS;}
    else{m=p;}
    
    dipl = compute_dipl(data.wfs(m).sX, data.tel.diam, hl, data.wfs(m).x/206265., data.wfs(m).y/206265., data.wfs(m).lgsH,xshift(m),yshift(m),obs=data.tel.obs);
    
    for(q=p; q<=nbWfs; q++) {       // boucle 2 sur les asos
      if(singleWFS){n = singleWFS;}
      else{n=q;}
      
      djql = compute_dipl(data.wfs(n).sX, data.tel.diam, hl, data.wfs(n).x/206265., data.wfs(n).y/206265., data.wfs(n).lgsH,xshift(n),yshift(n),obs=data.tel.obs);


      if(data.wfs(m).type!=0 && data.wfs(n).type!=0 ) {   // if both WFS are valid ones
       
        // subap size in the pupil plane
        sspSize1 = data.tel.diam/data.wfs(m).sX;
        sspSize2 = data.tel.diam/data.wfs(n).sX;
        
        tmp = Dij2covMat( dipl, djql, sspSize1, sspSize2,  cnh_hl, hl, l0h_hl,data.wfs(m).lgsH, data.wfs(n).lgsH,m,n,loworder=loworder,deriv=deriv);
        
        rp = tslrange(p);
        rq = tslrange(q);
        Cmod_hl(rp, rq) = tmp;

      }
    } // end loop WFS 2
  }   // end loop WFS 1

  return Cmod_hl;
}

func covMatAllLayers(nbWfs, cnh, alt, l0h,xshift,yshift, filtreTilt,loworder=,singleWFS=,deriv=,verb=)
/* DOCUMENT calcMatCov_explicit(nbWfs, L0, r0, alti, xaso, yaso, XPup, YPup, diamPup, thetaML)
   
   <nbWfs>      : scalar, number of wfss
   <L0>         : array 1D, scalar, log10(outer scale in meters), lenght = # of layers
   <r0>         : array 1D, lenght = # of layers
   <alti>       : array 1D, lenght = # of layers
   <filtreTilt> : tip-tilt is filtered on LGS when set to 1. Default=1.
   
   SEE ALSO:
   */
{
  
  if(!singleWFS){
    rr = sum(data.wfs.nssp); // total dim of the matrix
    Cmod = array(0.,rr,rr);
  }else{
    Cmod = array(0.,data.wfs(singleWFS).nssp,data.wfs(singleWFS).nssp);
  }

  // loop on turbulent layers
  nl = data.learn.nl;
  for(l=1; l<=nl; l++) {
    if(cnh(l) !=0){
      if(verb){
        write,format="\rGenerating the %d th layer...",l;
      }
      Cmod += covMat1layer(nbWfs, l,cnh(l) , alt(l), l0h(l),xshift,yshift,loworder=loworder,singleWFS=singleWFS,deriv=deriv);
    }
  }

  // Filtering of tiptilt according to the WFS-type : 0, 2 or 3 (no wfs, lgs, TTonly)

  if( is_void(filtreTilt ) ) filtreTilt = 1;
  if( filtreTilt==1 ) {
    Cmod = handle_tilt_from_wfstype(Cmod,ttr=data.learn.ttr);
  }else{
    if(!singleWFS)
      Cmod =  duplicateCovar(Cmod);
  }
  return Cmod;
  
}

func covMatModel(data, fitEstim,&grad,deriv=,loworder=,singleWFS=,verb=)
/* DOCUMENT 
     
     
   SEE ALSO:
 */
{
  local cnh,alt,l0h,tracking,xshift,yshift;
  // get all the variables that were encapsulated inside fitEstim
  unpackcoeffs, fitEstim, cnh, alt, l0h, tracking, xshift, yshift;
  
  Cfit = covMatAllLayers(data.nwfs, cnh, alt,l0h,xshift,yshift,filtreTilt,loworder=loworder,singleWFS=singleWFS);

   // ................... Tracking perturbation ........................
  if(!data.learn.ttr){
    Cfit = trackingMatCov(tracking, Cfit);
  }
  
  // 2-step case
  // this matrix is useful only when doing a "2-step Learn", for removing the ground layer
  if(data.learn.runLearnTwoSteps){
    P = *data.learn.transformationMatrix;
    Cfit = (P(,+) * Cfit(,+))(,+) *P(+,);
  }
  if(deriv){
    grad = array(0.,data.Nslopes,data.Nslopes,3*data.learn.nl);
    for(l=1;l<= data.learn.nl;l++){

      //bugs to be fixed: gradient about hl is not correct; TT on the L0 gradient
      // adds the gradient about the tracking, xshift,yshift (easy)
      
      //Cfit gradient about  L0
      grad(,,l) = handle_tilt_from_wfstype(covMat1layer(data.nwfs,l,cnh(l), alt(l), l0h(l),loworder=loworder,singleWFS=singleWFS,deriv=3));
      //Cfit gradient about alt
      grad(,,l+data.learn.nl) = handle_tilt_from_wfstype(covMat1layer(data.nwfs,l,cnh(l), alt(l), l0h(l),loworder=loworder,singleWFS=singleWFS,deriv=2));
      //Cfit gradient about cnh
      grad(,,l+2*data.learn.nl) = (covMat1layer(data.nwfs,l,cnh(l), alt(l), l0h(l),loworder=loworder,singleWFS=singleWFS,deriv=1));
      if( data.learn.diagonal==1) {
        takesDiag(grad(,,l)) = 0.;
      }
    }
  }
  // ................... Noise options (diagonal or not)  .......................
  if( data.learn.diagonal==1) {
    takesDiag(Cfit) = 0.;
  }

  return Cfit;
}

func trackingMatCov(tracking,Cfit)
/* DOCUMENT trackingMatCov(tracking, fullMat)
     Creates a checkerboard on top of the covariance matrix, with
     <x^2>=tracking(1) on all the x/x covariance areas,
     <y^2>=tracking(2) on all the y/y covariance areas,
     <xy> =tracking(3) on all the x/y and y/x covariance areas.
   SEE ALSO:
 */
{
  fullMat = Cfit;
  for(p=1;p<=data.nwfs;p++) {
    for(q=p;q<=data.nwfs;q++) {
      if( data.wfs(p).type!=2 && data.wfs(q).type!=2 ) {
        // ranges for X-slopes
        xrp = tslindex(p):(tslindex(p)+data.wfs(p).nssp/2-1);
        xrq = tslindex(q):(tslindex(q)+data.wfs(q).nssp/2-1);
        // ranges for Y-slopes
        yrp = tslindex(p)+data.wfs(p).nssp/2:tslindex(-p);
        yrq = tslindex(q)+data.wfs(q).nssp/2:tslindex(-q);
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
  // 2*gamma_R(11./6)*2^(-5./6)*pi^(-8./3)*(24*gamma_R(6./5)/5.)^(5./6);
  k1 = 0.1716613621245709486;
  dprf0 = (2*pi/L0)*r;
  // k2 is the value for gamma_R(5./6)*2^(-1./6),
  // but is now unused
  // k2 = 1.0056349179985892838;

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
  Ga  =  2.01126983599717856777;       // Gamma(a) / (1/2)^a
  Gma = -3.74878707653729348337;       // Gamma(-a) * (1/2.)^a
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
  // k2 is the value for
  // gamma_R(5./6)*2^(-1./6)
  k2 = 1.00563491799858928388289314170833;
  k3 = 1.25331413731550012081;   //  sqrt(pi/2)
  a1 = 0.22222222222222222222;   //  2/9
  a2 = -0.08641975308641974829;  //  -7/89
  a3 = 0.08001828989483310284;   // 175/2187
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
/*
  ____               _ _          
 / ___|_ __ __ _  __| (_) ___ _ __ | |_ 
| |  _| '__/ _` |/ _` | |/ _ \ '_ \| __|
| |_| | | | (_| | (_| | |  __/ | | | |_ 
 \____|_|  \__,_|\__,_|_|\___|_| |_|\__|                                      
*/

func gradDPHI(x,y,cnh,L0,p,q,deriv=)
{
   if(cnh==0 | L0==0) return 0*x;
  mod = "Von-Karman";
  if(L0 >=1e4) mod = "Kolmogorov";
    
  if(deriv == 1){
    return DPHI(x,y,L0)/abs(cnh);
  }
  else if(deriv == 2){
    if(p==q) return 0*x;
    return grad_Dphionr(x,y,L0,model=mod) * grad_ronhl(x,y,p,q);
  }else if(deriv == 3){
    if(mod == "Von-Karman")
      return grad_DphionL0(x,y,L0);
    else return 0*x;
  }
}
func grad_DphionL0(x,y,L0)
{
  r = abs(x,y);
  
  //4*pi*pi*k1 =  6.77692
  //5*k1*k2/3 = 0.287714
  grad =  0.287714*L0^(2/3.) - 6.77692*r*r*(L0^(-4/3.))*Kalpha(2*pi*r/L0,-1/6.);

  return grad;
}

func grad_ronhl(x,y,p,q)
{

  grad = 0*x;
  nx = dimsof(x)(0);
  grad_x = array(0.,nx);
  grad_y = array(0.,nx);

  rix = riy = rjx = rjy = 0;
  if(data.wfs(p).lgsH){
    ri = compute_ri(data.wfs(p).sX, data.tel.diam, data.tel.obs)/data.wfs(p).lgsH;
    rix = ri(,1);
    riy = ri(,2);
  }
  if(data.wfs(q).lgsH){
    rj = compute_ri(data.wfs(q).sX, data.tel.diam, data.tel.obs)/data.wfs(q).lgsH;
    rjx = rj(,1);
    rjy = rj(,2);
  }

  grad_x += (data.wfs(p).x - data.wfs(q).x)/206265. + rjx - rix;
  grad_y += (data.wfs(p).y - data.wfs(q).y)/206265. + rjy  - riy;

  for(i=1;i<=nx;i++){
    grad(i,) = grad_x(i);
    grad(,i) = grad_y(i);
  }

  return grad;
}

func grad_Dphionr(x,y,L0,model=)
{
  if(is_void(model)) model = "Von-Karman";

  r = abs(x,y)+1e-15;
  
  if(model == "Von-Karman"){
    //4*pi*pi*k1 =  6.77692
    grad =  6.77692*r*(L0^(-1./3))*Kalpha(2*pi*r/L0,-1/6.);
  }else if(model == "Kolmogorov"){
    grad = 5*DPHI(x,y,L0,model = "Kolmogorov")/(3*r);
  }
  return grad;
}

func Kalpha(x,alpha,k=)
/* DOCUMENT f(x) = x^(alpha)*Kalpha(x)
   
     Computation of the function
     f(x) = x^(alpha)*K_{alpha}(x)
     using a series for the esimation of K_{alpha}, taken from Rod Conan thesis :
     K_alpha(x)=1/2 \sum_{n=0}^\infty \frac{(-1)^n}{n!}
     (\Gamma(-n-alpha) (x/2)^{2n+alpha} + \Gamma(-n+a) (x/2)^{2n-alpha}) ,
     

     Setting x22 = (x/2)^2, setting uda = (1/2)^alpha, and multiplying by x^alpha,
     this becomes :
     x^a * Ka(x) = 0.5 $ -1^n / n! [ G(-n-a).uda x22^(n+a) + G(-n+a)/uda x22^n ]
     Then we use the following recurrence formulae on the following quantities :
     G(-(n+1)-alpha) = G(-n-alpha) / -alpha-n-1
     G(-(n+1)+alpha) = G(-n+alpha) /  alpha-n-1
     (n+1)! = n! * (n+1)
     x22^(n+1) = x22^n * x22
     and at each iteration on n, one will use the values already computed at step (n-1).

     The first term of the series has also been skipped, as it
     vanishes with another term in the expression of Dphi.
     
   SEE ALSO:
 */
{
  if( is_void(k) ) k=10;
  fn = 1.;                             // initialisation factorielle 0!=1
  w = where(x !=0);
  x2a = x;
  x2a(w) = (x(w))^(2.*alpha);
  x22 = x*x/4.;                        //  (x/2)^2
  x2n = 0.5;                           // init (1/2) * x^0
  
  Ga  =  gamma_R(alpha)/.5^(alpha);       // Gamma(a) / (1/2)^a
  Gma = gamma_R(-alpha)*.5^(alpha);      // Gamma(-a) * (1/2.)^a
  
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
      Gma /= -alpha-n-1; // gamma function
      Ga /= alpha-n-1;   // idem
      x2n *= x22;    // x^n
    }
  }
  return s;
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
  
  MatFiltZer = identite(data.Nslopes);
  if(!filt_ngs) filt_ngs=0;
  if(!filt_lgs) filt_lgs=0;
  if(!filt_ttgs) filt_ttgs=0;
  
  if(modefilt!=0){
    
    for(i=1; i<=data.nwfs; i++) {
      
      if((data.wfs(i).type==1 && filt_ngs==1) || (data.wfs(i).type==2 && filt_lgs==1)  || (data.wfs(i).type==3 && filt_ttgs==1) ) {
        mrz = *data.wfs(i).mrz;
        mrz(modefilt+1:,) = 0;

        zmi = (*data.wfs(i).zmi)(,1:modefilt);
        mrz = invgen(zmi,0, t=1);
        nslopes_i = data.wfs(i).nssp;
        MatFiltZer(slrange(i),slrange(i)) = identite(nslopes_i) - (zmi)(,+) * mrz(+,);
      }
      // si ASO desactive on vire tout
      if(data.wfs(i).type==0 || (data.wfs(i).type==1 && filt_ngs==-1) || (data.wfs(i).type==2 && filt_lgs==-1)  || (data.wfs(i).type==3 && filt_ttgs==-1))
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

    MatKeepZer = identite(data.Nslopes) - Filt_Zer(modefilt,filt_ngs=calc_ngs,filt_lgs=calc_lgs,filt_ttgs=calc_ttgs);
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
  MatCompGL = array(0.,data.Nslopes,data.Nslopes);
  if(!avg_ngs) avg_ngs=0;
  if(!avg_lgs) avg_lgs=0;
  if(!avg_ttgs) avg_ttgs=0;
  
  activeWfs = (data.wfs.type==1 & avg_ngs==1) | (data.wfs.type==2 & avg_lgs==1) | (data.wfs.type==3 & avg_ttgs==1);
  
  nbActiveWfs = 0;
  for(i=1;i<=data.nwfs;i++) {
    if( activeWfs(i) ) {
      nbActiveWfs++;
      for(j=1;j<=data.nwfs;j++) {
        if( activeWfs(j) ) {
          MatCompGL(slrange(i),slrange(j)) = MatCompGL(slrange(j),slrange(i)) = identite(data.wfs(i).nssp);
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
    P =  identite(data.Nslopes) - (MatCompGL_TT(,+)*MatKeepTT(,+) + MatCompGL_HO(,+)*MatKeepHO(,+));
  }else{
    //We only average slopes that are all tip-tilt removed
    MatKeepHO = Filt_Zer(modefilt,filt_ngs=1,filt_lgs=1,filt_ttgs=-1);
    MatCompGL_HO = Calc_CommonMat(avg_ngs=1,avg_lgs=1,avg_ttgs=0);
    P =  identite(data.Nslopes) -  MatCompGL_HO(,+)*MatKeepHO(,+);
  }
  
  return P;

}





func duplicateCovar( caa_input )
{

  for(i=1; i<=data.nwfs; i++) {
    for(j=i; j<=data.nwfs; j++) {
      
      ri = tslrange(i);
      rj = tslrange(j);
      
      if( data.wfs(i).type==0 || data.wfs(j).type==0 ) {
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
  
  for(i=1; i<=data.nwfs; i++) {
    for(j=i; j<=data.nwfs; j++) {
      
      ri = tslrange(i);
      rj = tslrange(j);
      
      tmp = caa_full(ri,rj);
      if(!ttr){
        if(data.wfs(i).type==1) {
          tmpi = tmp;
        } else if(data.wfs(i).type==2) {
          tmpi = (*data.wfs(i).filterouttiltMatrix)(,+)*tmp(+,);
        } else if(data.wfs(i).type==3) {
          tmpi = (*data.wfs(i).keeponlytiltMatrix)(,+)*tmp(+,);
        } else {
          tmpi = tmp * 0.0;
        }
        if(data.wfs(j).type==1) {
          caa_full(ri,rj) = tmpi;
        } else if(data.wfs(j).type==2) {
          caa_full(ri,rj) = tmpi(,+) * (*data.wfs(j).filterouttiltMatrix)(+,);
        } else if(data.wfs(j).type==3) {
          caa_full(ri,rj) = tmpi(,+) * (*data.wfs(j).keeponlytiltMatrix)(+,);
        } else {
          caa_full(ri,rj) = tmpi * 0.0;
        }
      }else{
        tmpi = (*data.wfs(i).filterouttiltMatrix)(,+)*tmp(+,);
        caa_full(ri,rj) = tmpi(,+) * (*data.wfs(j).filterouttiltMatrix)(+,);
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










/*
  if(useRawMat){
      if(verb){
        write,"\rDerivating the tomographic residual covariance matrix...";
      }
      //computes the raw  covariance matrix from the tomographic distance
      Cdd = (*data.slopes_dis)(,+) * (*data.slopes_dis)(,+)/data.nframes;
      Cee = computeCovSlopesError(Cdd,*data.rtc.R);
      //unbiases from the TS noise 
      Cee -=  Cnn(slrange(data.its),slrange(data.its));
      //manages the time-filtered off-axis noise contribution
      noisefact = filteringNoiseFactor(data.rtc.gain,data.rtc.delay,data.rtc.obsmode);
      if(verb){
        write,format="\rNoise filtering factor: %g\n",noisefact;
      }
      Cee += noisefact * computeCovSlopesError(Cnn,*data.rtc.R);
      //manages the aliasing
      Caa = computesCaaMatrix(data.nwfs,data.turbu.cnh,data.turbu.altitude,data.turbu.l0h,data.wfs.x,data.wfs.y,verb=verb);
      Cee += (*data.rtc.R)(,+) * ((Caa(norange(data.its),norange(data.its)))(,+)*(*data.rtc.R)(,+))(+,) - computeCovSlopesError(Caa,*data.rtc.R);
    }

*/







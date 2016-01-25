func simulAOsystem(void,verb=)
/* DOCUMENT

 */
{

  //constants
  arc2nm  = 0.5*tel.diam*1e9/radian2arcsec;
  nm2rad  = 2*pi*1e-9/cam.lambda;
  
  //////////////////////////////////
  // .... DEFINING GEOMETRY .... //
  ////////////////////////////////

  N     = tel.nPixels;
  n     = 256;
  x     = (indgen(N) -N/2)(,-:1:N);
  x    *= tel.pixSize/(tel.diam/2.);
  y     = transpose(x);
  tmp   = polarFromCartesian(x,y);
  rho   = tmp(,,1);
  theta = tmp(,,2);

  //Generating and saving zernike modes
  Zi = readfits("fitsFiles/zernikeModes.fits",err=1);
  if(is_void(Zi)){
    Zi = array(0.0,N,N,35);
    for(i=1;i<=35;i++){
      Zi(,,i) = computeZernikePolynomials(rho,theta,i);
    }
    writefits,"fitsFiles/zernikeModes.fits",Zi;
  }

  P = circularPupFunction(x,y,tel.obs);
  
  /////////////////////////////////
  // .... LOADING TELEMTRY .... //
  ///////////////////////////////

  // ..... Loading matrices
  mc    = *rtc.mc; // in volts/arcsec
  if(rtc.obsMode == "MOAO"){
    R     = *rtc.R; // to be applied on arcsec
    //command matrix
    mct = mc(,+) * R(+,); // in volts/arcsec
  }
  //Zernike reconstruction matrix
  mi   = *rtc.mi // in arcsec/volts
  mrz  = *sys.slopesToZernikeMatrix;
  v2z  = mrz(,+) * mi(+,) *  arc2nm * nm2rad
  nzer = dimsof(mrz)(2);
  
  // ..... Loading slopes
  slopes = *rtc.slopes_dis; // in arcsec
  Son   = slopes(slrange(rtc.its),);
  Soff  = slopes(norange(rtc.its),);
  //Son  -= Son(,avg);
  //Soff -= Soff(,avg);

  // ..... Tomography
  if(rtc.obsMode == "MOAO"){
    Utomo    = mct(,+) * Soff(+,); // in volts
  }else{
    Uon      = mc(,+) * Son(+,);
  }
  
  // Atmospheric phase
  a_atm   = mrz(,+) * Son(+,);
  a_atm   *= arc2nm * nm2rad; // in radian
  a_dm    = 0*a_atm;
  //allocation in memory
  a_res   = array(0.,nzer,niter);
  phi_res = psf_cam = array(0.,N,N);
  SR_inst = SR_avg = array(0.,niter);
  
  ///////////////////////////////
  // .... LAW OF COMMAND .... //
  //////////////////////////////

  // ..... Loading loop parameters
  niter    = 128;//dimsof(slopes)(0);
  nActu    = dm.nValidActu;
  loopGain = rtc.loopGain;
  dt       = rtc.frameDelay%1;
  //instantation of voltages parameters
  U =  V = array(0.,nActu,niter);
  // offset volatges
  V0 = 0;
  suff = readFitsKey(pathdata,"STATIC");
  V0 = restorefits("offsetDM",suff);//offset in volts
  V0 = V0(1:-2);

  Nexp = int(cam.exposureTime*rtc.Fe);
  n = 0;
  nzer = dimsof(*sys.slopesToZernikeMatrix)(2);
  ares = array(0.,nzer);
  
  if(rtc.obsMode == "MOAO"){
    for(k=3;k<=niter;k++){
      //index for exposure time
      n++; a_res *=0;
      //defining new voltages
      U(,k)     = (1 - loopGain) * U(,k-1) + loopGain*Utomo(,k);
      V(,k)     = U(,k) + V0;

      //defining parallel phase
      a_dm(,k)  = v2z(,+) * V(+,k); // in radian

      //delayed command
      a_delay = a_dm(,k-2)*dt + (1-dt)*a_dm(,k-1);
        
      //residue
      a_res += a_atm(,k) -  a_delay;// perp phase ??

      
      //residual phase
      for(i=1;i<=nzer;i++){
        phi_res += a_res(i) * Zi(,,i); // in radians
      }

      //defining instantaneous psf
      /*E = P*exp(1i*phi_res);
      //error;
      psf_inst  = roll(abs(fft(roll(E)))^2);
      psf_inst /= sum(psf_inst);
      psf_inst /= tel.airyPeak;
      
      psf_cam  += psf_inst;
      psf_cam  *= (k-1.)/(k);

      SR_inst(k) = max(psf_inst);
      SR_avg(k)  = max(psf_cam);
      */
      
      //display
      if(disp){
        window,10;clr;
        pli,psf_cam,-N*tel.pixSize/2,-N*tel.pixSize/2,N*tel.pixSize/2,N*tel.pixSize/2;
        limits,-n*tel.pixSize,n*tel.pixSize;
        range, -n*tel.pixSize,n*tel.pixSize;
        xytitles,"arcsec","arcsec";
      }
      write,format="Job done:%.3g%s\r",100.*k/niter,"%";
    }
    
  }else if(rtc.obsMode == "SCAO"){
    for(k=3;k<=niter;k++){
      //defining new voltages
      U(,k)  = U(,k-1) + loopGain*Uon(,k);

      //defining parallel phase
      a_dm(,k)  = v2z(,+) * U(+,k);

      //delayed command
      a_delay = a_dm(,k-2)*dt + (1-dt)*a_dm(,k-1)
      //residue
      a_res = a_atm(,k) - a_delay;//add the perp phase
      
      //residual phase
      for(i=1;i<=nzer;i++){
        phi_res += a_res(i) * Zi(,,i);
      }

      //defining instantaneous psf
      E = P*exp(1i*phi_res);
      psf_inst  = roll(abs(fft(roll(E)))^2);
      psf_inst /= sum(psf_inst);
      psf_inst /= tel.airyPeak;

      psf_cam  += psf_inst;
      psf_cam  *= (k-2.)/(k-1.);

      SR_inst(k) = max(psf_inst);
      SR_avg(k)  = max(psf_cam);
      
      //display
      if(disp){
        window,10;clr;
        pli,psf_cam,-N*tel.pixSize/2,-N*tel.pixSize/2,N*tel.pixSize/2,N*tel.pixSize/2;
        limits,-n*tel.pixSize,n*tel.pixSize;
        range, -n*tel.pixSize,n*tel.pixSize;
        xytitles,"arcsec","arcsec";
      }
      write,format="Job done:%.3g%s\r",100.*k/niter,"%";
        
    }

  }


  window,10;clr;
  pli,psf_cam,-N*tel.pixSize/2,-N*tel.pixSize/2,N*tel.pixSize/2,N*tel.pixSize/2;
  limits,-n*tel.pixSize,n*tel.pixSize;
  range, -n*tel.pixSize,n*tel.pixSize;
  xytitles,"arcsec","arcsec";

  return [SR_inst,SR_avg];       
}

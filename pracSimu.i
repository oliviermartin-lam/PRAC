include, "pracMain.i"; 



func pracSimu(verb=,disp=)
/*DOCUMENT p = pracSimu(verb=1,disp=1);

 
 */
{


  include, "pracConfig.i",1;
  geometry = "square";

  
  /////////////////////
  //.... Initializing the data struct
  //////////////////////////////////////////////////////////
  
  include, "simulStructConfig.i",1;
  defineSimulStructs,verb=verb;

  /////////////////////
  // .... perfect telescope OTF
  //////////////////////////////////////////////////////////
  
  OTF_tel = OTF_telescope(tel.diam,tel.obs,tel.nPixels,tel.pixSize);
    
  /////////////////////
  // .... OTF from fitting error
  //////////////////////////////////////////////////////////
   
  OTF_fit = computeOTFfitting(geometry,verb=verb);
     
  /////////////////////
  // .... OTF from response delay time of the system
  //////////////////////////////////////////////////////////
  
  OTF_bw   = computeOTFbandwidth(geometry,rtc.obsMode,verb=verb);

  /////////////////////
  // .... OTF from tomographic residue
  //////////////////////////////////////////////////////////////////////////////

 
  // Computing the tomographic residue cov Matrix.
  Cee = computeCovSlopesError(*covMatrix.atm,*rtc.R);

  //computes the covariance matrix of voltages in volts^2
  Cvv = propagateError2Dm(Cee, *rtc.mc );

  // computation of influence function and defining the pupil
  N = tel.nPixels;
   
  //Getting the map
  Cvvmap = getMap(Cvv,1);
  ncov = dimsof(Cvvmap)(0);// nber of elements on each side of the center of Cvvmap
  ncmap = N/2+1;
  nelem = (ncov-1)/2;

  //defining the new map in the real domain
  start = int(ncmap - tel.diam/tel.pixSize);
  stop  = int(ncmap + tel.diam/tel.pixSize);
  step  = tel.pitchSize;
  rr    = start:stop:step;

  cv = array(0.0,N,N);
  cv(rr,rr) = Cvvmap; 

  //TF of the influence function autocorrelation
  x  = (indgen(N) -N/2-1)(,-:1:N);
  x *= tel.pixSize/(tel.diam/2.);
  y  = transpose(x);
  fi = funcInflu(x,y,dm.x0);
  A  = abs(fft(fi))^2; // in m^-2
    
  //computing the correlation function
  Sp = numberof(fi)*tel.pitchSize^2;
  C  = roll(fft( A * fft(roll(cv)) ).re / Sp); // in m^2
    
  //getting the phase structure function in rd^2
  dphi  = max(C) - C; 
  dphi  = roll(dphi); 
  dphi *= (2*pi/cam.lambda)^2.;

  //computinh the OTF
  OTF_tomo  = exp(-0.5*dphi);
    
       
  /////////////////////
  // .... OTF at the TS location
  //////////////////////////////////////////////////////////
  
  OTF_res =  OTF_tel * OTF_fit * OTF_bw * OTF_tomo ;

       
  //telescope included into OTF_stats

  PSF_res = roll( fft(OTF_res).re );

  
  /////////////////////
  // .... PERFORMANCE
  //////////////////////////////////////////////////////////

  FWHM_res = getPsfFwhm(PSF_res,tel.pixSize,fit=2);
  SR_res   = sum(OTF_res);
  SR_fit   = sum(OTF_fit  * OTF_tel);
  SR_bw    = sum(OTF_bw   * OTF_tel);
  SR_tomo  = sum(OTF_tomo * OTF_tel);
  
  /////////////////////
  // .... Display and verbose
  //////////////////////////////////////////////////////////

  if(disp){
    l = tel.nPixels * tel.pixSize;
    window,1; clr;pli, PSF_res,-l/2,-l/2,l/2,l/2;
    limits,-1,1;range,-1,1;
    pltitle,"Estimated PSF with SR = " + var2str(arrondi(100*SR_res,1))+"%";
    xytitles,"Arcsec","Arcsec";
  }
  
  if(verb){
    write,format="SR reconstructed    = %.3g%s\n", 100*SR_res,"%";
    write,format="FWHM reconstructed  = %.4g%s\n", 1e3*FWHM_res," mas";
    write,format="SR fit              = %.3g%s\n", 100*SR_fit,"%";
    write,format="SR bw               = %.3g%s\n", 100*SR_bw,"%";
    write,format="SR tomo+alias+noise = %.3g%s\n", 100*SR_tomo,"%";    
  }

  return PSF_res;
}

func makeOTF(Cvv,&dphi,averageMode=,verb=)
/* DOCUMENT

 */
{
  tic,10;
  
  if(verb)
    write,format="%s ... ","Interpolating Dphi map";

  if(is_void(averageMode)) averageMode = "Vii";
  
  // computation of influence function and defining the pupil
  N = tel.nPixels;
  P = *tel.pupil;

  //autocorrelation of the pupil expressed in pupil radius^-1
  P2         = enlargeSupport(P,tel.nTimes);
  fftpup     = fft(roll(P2));
  conjpupfft = conj(fftpup);
  G          = fft(fftpup*conjpupfft).re;
  //Interpolating
  //conjpupfft = interpolateOTF(conjpupfft.re,N) + 1i*interpolateOTF(conjpupfft.im,N);
  //G          = interpolateOTF(G,N);
  //defining the inverse
  den             = 0*G;
  msk             = G/max(G) > 1e-5;
  den(where(msk)) = 1./G(where(msk));
  
  //defining the modes
  fi = *dm.modes;

  if(averageMode == "Uij"){

    /*
       _   _ ___    _ 
      | | | |_ _|  | |
      | | | || |_  | |
      | |_| || | |_| |
       \___/|___\___/ 
                
  */
    
    tmp = 0*P2;
    for(i=1;i<=dm.nValidActu;i++){
      for(j=i;j<=dm.nValidActu;j++){
        //modes
        Mi   = fi(,,i);
        Mj   = fi(,,j);
        //Enlarging support
        Mi   = roll(enlargeSupport(Mi,tel.nTimes));
        Mj   = roll(enlargeSupport(Mj,tel.nTimes));
        //Uij computation
        t1   = (fft(Mi*Mj)*conjpupfft).re;
        t2   = (fft(Mi)*conj(fft(Mj))).re;
        Uij  = fft(2*(t1 - t2),[-1,-1]).re;
        //summing modes
        tmp +=  ((i!=j) + 1.)*Cvv(i,j)*Uij;
      }
    }
    //multiplying by a mask
    dphi = den * tmp * (2*pi/cam.lambda)^2.;
    //to get the SF in rd^2
    OTF  = exp(-0.5*dphi);
    
  }else if(averageMode == "Vii"){


    /*
      __     _____ ___ 
      \ \   / /_ _|_ _|
       \ \ / / | | | | 
        \ V /  | | | | 
         \_/  |___|___|
                 
    */

   
    //Diagonalizing the Cvv matrix
    l    = SVdec(Cvv,Bt,B); // in volts^2
    nll  = dm.nValidActu;
    cond = 100.;
    nll  = numberof(where(l>max(l)/cond))
    //Cvv = B*l(,-))(+,)*Bt(,+)
    M = fi(,,+)*B(,+); 
  
    //loop on actuators
    tmp = Mi = 0.*P2;
    for(i=1;i<=nll;i++){
      Mi   = M(,,i)*P; //must be in volt^2 * meter^2
      Mi   = roll(enlargeSupport(Mi,tel.nTimes));
      //Vii computation  in m^-2.volts^-2
      Vii  = (fft(Mi*Mi)*conjpupfft).re - abs(fft(Mi))^2;
      //summing modes into dm basis
      tmp += l(i) * Vii; //must be in meters^-2
    }

    dphi = den * fft(2*tmp,[-1,-1]).re * (2*pi/cam.lambda)^2. ; 
    OTF  = exp(-0.5*dphi);

  }else if(averageMode == "intersample"){


    /*
      ___ _   _ _____ _____ ____  ____    _    __  __ ____  _     _____ 
     |_ _| \ | |_   _| ____|  _ \/ ___|  / \  |  \/  |  _ \| |   | ____|
      | ||  \| | | | |  _| | |_) \___ \ / _ \ | |\/| | |_) | |   |  _|  
      | || |\  | | | | |___|  _ < ___) / ___ \| |  | |  __/| |___| |___ 
     |___|_| \_| |_| |_____|_| \_\____/_/   \_\_|  |_|_|   |_____|_____|
                                                                   
    */
    
    //Getting the map
    Cvvmap = getMap(Cvv,1);
    ncov   = dimsof(Cvvmap)(0);// nber of elements on each side of the center of Cvvmap
    ncmap  = N/2+1;
    nelem  = (ncov-1)/2;

    //defining the new map in the real domain
    start = int(ncmap - tel.diam/tel.pixSize);
    stop  = int(ncmap + tel.diam/tel.pixSize);
    step  = tel.pitchSize;
    rr    = start:stop:step;
    
    cv        = array(0.0,N,N);
    cv(rr,rr) = Cvvmap; 

    //TF of the influence function autocorrelation
    tmp = pupilGeometry(N);
    fi = roll(funcInflu(tmp(,,1),tmp(,,2),dm.x0));
    A  = abs(fft(fi))^2; // in m^-2
    
    //computing the correlation function
    Sp = numberof(fi)*tel.pitchSize^2;
    //Sp = numberof(fi)/double(N/dm.nActu)^2;
    C  = fft( A * fft(roll(cv)) ).re / Sp; // in m^2
    
    //getting the phase structure function in rd^2
    dphi  = max(C) - C; 
    dphi *= (2*pi/cam.lambda)^2.;

    //computing the OTF
    OTF  = exp(-0.5*dphi);

  }

  if(verb)
    write,format="done in %.3g s\n",tac(10);

  //error;
  //OTF = roll(enlargeSupport(roll(OTF*msk),tel.nTimes));
    
  return OTF*msk;
}

func makeOTF(Cvv,&dphi,averageMode=,verb=)
/* DOCUMENT

 */
{
  tic,3;
  
  if(verb)
    write,format="%s ... ","Interpolating Dphi map";

  if(is_void(averageMode)) averageMode = "Vii";
  
  // computation of influence function and defining the pupil
  N = tel.nPixels;
  P = circularPupFunction(tel.diam,tel.obs,tel.nPixels,tel.pixSize) ;

  //autocorrelation of the pupil expressed in pupil radius^-1
  fftpup = fft(P);
  conjpupfft = conjugate(fftpup);
  G = fft(fftpup*conjpupfft).re;
  //defining the inverse
  den = 0*G;
  msk = G/max(G) > 1e-5;
  den(where(msk)) = 1./G(where(msk));
  
  //defining the modes
  if(averageMode != "intersample"){
    fi = readfits("fitsFiles/dm_modes.fits",err=1); // influence functions in meters
    if(is_void(fi)){
      fi = array(0.,N,N,dm.nValidActu);
      xi = (*dm.csX);
      yi = (*dm.csY);
      x     = (indgen(N) -N/2)(,-:1:N);
      x    *= tel.pixSize/(tel.diam/2.);
      y     = transpose(x);
      for(i=1;i<=dm.nValidActu;i++){
        fi(,,i) = funcInflu(x-xi(i),y-yi(i),dm.x0);
      }
      fi *= double(P(,,-));
      writefits,"fitsFiles/dm_modes.fits",fi;
    }
  }else{
    x     = (indgen(N) -N/2)(,-:1:N);
    x    *= tel.pixSize/(tel.diam/2.);
    y     = transpose(x);
    fi = funcInflu(x,y,dm.x0);
  }

  if(averageMode == "Uij"){
    tmp = 0*P;
    for(i=1;i<=dm.nValidActu;i++){
      for(j=i;j<=dm.nValidActu;j++){
        //modes
        Mi   = fi(,,i);
        Mj   = fi(,,j);
        //Uij computation
        t1   = (fft(Mi*Mj)*conjpupfft).re;
        t2   = (fft(Mi)*conjugate(fft(Mj))).re;
        Uij  = fft(t1 - t2);
        //summing modes
        tmp +=  ((i!=j) + 1.)*Cvv(i,j)*Uij;
      }
    }
    //multiplynig by a mask
    dphi = den * tmp.re * (2*pi/cam.lambda)^2.;
    //to get the SF in rd^2
    OTF  = exp(-0.5*dphi);

  }else if(averageMode == "Vii"){
    
    //Diagonalizing the Cvv matrix
    l = SVdec(Cvv,Bt,B); // in volts^2
    //Cvv = B*l(,-))(+,)*Bt(,+)
    M = fi(,,+)*B(,+); 
  
    //loop on actuators
    tmp = Mi = 0.*P;
    for(i=1;i<=dm.nValidActu;i++){
      Mi   =  M(,,i); //must be in volt^2 * meter^2
      //Vii computation  in m^-2.volts^-2
      Vii  = (fft(Mi*Mi)*conjpupfft).re - abs(fft(Mi))^2;
      //summing modes into dm basis
      tmp += l(i) * Vii; //must be in meters^-2
    }

    dphi = den * fft(tmp).re * (2*pi/cam.lambda)^2. ; 
    OTF  = exp(-0.5*dphi);

  }else if(averageMode == "intersample"){

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
    A = abs(fft(fi))^2; // in m^-2
    
    //computing the correlation function
    Sp = numberof(fi)*tel.pitchSize^2;
    C  = roll(fft( A * fft(roll(cv)) ).re / Sp); // in m^2
    
    //getting the phase structure function in rd^2
    dphi  = max(C) - C; 
    dphi  = roll(dphi); 
    dphi *= (2*pi/cam.lambda)^2.;

    //computinh the OTF
    OTF  = exp(-0.5*dphi);
  }

  if(verb)
    write,format="done in %.3g s\n",tac(3);

  
  return msk*OTF;
}

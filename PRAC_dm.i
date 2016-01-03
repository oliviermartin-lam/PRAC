/*
 ___       _                      _   _                               _ 
|_ _|_ __ | |_ ___ _ __ __ _  ___| |_(_) ___  _ __     __ _ _ __   __| |
 | || '_ \| __/ _ \ '__/ _` |/ __| __| |/ _ \| '_ \   / _` | '_ \ / _` |
 | || | | | ||  __/ | | (_| | (__| |_| | (_) | | | | | (_| | | | | (_| |
|___|_| |_|\__\___|_|  \__,_|\___|\__|_|\___/|_| |_|  \__,_|_| |_|\__,_|
                                                                        
  ____                                          _ 
 / ___|___  _ __ ___  _ __ ___   __ _ _ __   __| |
| |   / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` |
| |__| (_) | | | | | | | | | | | (_| | | | | (_| |
 \____\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|
                                                  
 __  __       _        _               
|  \/  | __ _| |_ _ __(_) ___ ___  ___ 
| |\/| |/ _` | __| '__| |/ __/ _ \/ __|
| |  | | (_| | |_| |  | | (_|  __/\__ \
|_|  |_|\__,_|\__|_|  |_|\___\___||___/
                                       
*/

func intermat( dm )
/* DOCUMENT mia = intermat( dm )
   
     Returns the interaction matrix, in arcseconds/volt.
     
   SEE ALSO:
 */
{

  // truth sensor number ...
  its = data.its;
  
  // Creation of subaperture points of the TS
  nssp = data.wfs(its).sX;
  xy =  compute_ri(nssp, 2.0, data.tel.obs,0,0);

  // matrix of distances (in x and y) from subaperture centre to DM
  dx = xy(,1) - (*data.dm.csX)(-,);
  dy = xy(,2) - (*data.dm.csY)(-,);

  // subap size
  ssize = 2./nssp;
  
  // creation of MI
  nsubap = int(data.wfs(its).nssp/2);
  mia = array(0.0, nsubap*2, data.dm.nactu);

  // here <mia> has no unit, it is just a difference between the edges of a subap
  // The result is expressed in meters.
  mia(1:nsubap,)  = difffuncInflu(dx, dy, ssize, data.dm.x0);
  mia(nsubap+1:,) = difffuncInflu(dy, dx, ssize, data.dm.x0);

  // here we convert mia into an angle (in rad) per volt
  // We need to divide the OPD by the size of subap to get an angle (radians)
  mia /= data.tel.diam/nssp;
  //we switch from radians to arcseconds
  mia *= 206265.;
  
  return mia;
  
}

func difffuncInflu(x,y,ssize,x0)
/* DOCUMENT d_funcInflu(x,y,ssize,x0)

   Returns the difference at the edge of the subaperture. The
   variables x, y, ssize and x0 are all in the same units. The
   returned value does not depend on the chosen unit, as it is just a
   difference between the 2 edges.

   This difference is an optical path difference, expressed in METERS.
   
   SEE ALSO:
 */
{
  ss = ssize / 2;
  return funcInflu(x + ss,y,x0) - funcInflu(x - ss,y,x0);
}

func funcInflu(x,y,x0)
/* DOCUMENT opd_metres = funcInflu(x,y,x0)

   The arguments <x>, <y> and <x0> must all be in the same units.
   
   In the particular case where <x0> is set to x0=dm.x0, <x> and <y>
   must be expressed with no unit, they are variables defined over the
   unitary circle of R=1.

   In any case, the function returns the influence function (=OPD,
   optical path difference) of the actuators expressed in METERS.
     
   SEE ALSO:
 */
{
  // allez on va dire que ce sont des metres !
  return 1e-6 * exp( -(x*x+y*y)/(2*x0*x0) );
}


func computeCommandMat( mia, threshold=, nmf=, condi=, win=,disp= )
{
 
  if( nmf==[] )
    nmf = -1;
  if( threshold==[] )
    threshold = 0.37;
    
  s = dimsof(mia);
  i = s(2);  // number of slopes
  j = s(3);  // number of actuators

 
  // SVD
  l = SVdec(mia, U, VT);
    
  // check eigen values ..........
  // 
  check_eigen, l, nmf, "command matrix eigenvals", disp=disp, thres=threshold, condi=sqrt(condi), win=win;   // nmf is an output !!
    
  // inverting eigen values
  l1 = invert_eigen( l, nmf );
    
  // computing the pseudo-inverse
  // mca = (U*l1(-,))(,+) * VT(+,);  // output has same shape as input 
  mca = VT(+,) * (U*l1(-,))(,+);     // output has trasposed shape wrt input

  return mca;
}
  
func invgen(mat , nfilt, t=)
/* DOCUMENT invgen( mat, nfilt, t=1 )
	mat   : matrix to be inverted
	nfilt : number of modes to be filtered
        t=1   : output has transposed shape wrt input
        t=0 (default) : output has same shape as input
*/
{
    s = SVdec(mat,u,vt);
    if( nfilt>0 ) {
        s1 = s;
        s1(1-nfilt:0) = 1.0;
        s1 = 1./s1;
        s1(1-nfilt:0) = 0.0;
    } else {
        s1 = 1.0 / s;
    }
    if( t==1 )
      m1 =  vt(+,) * (u*s1(-,))(,+); // output has trasposed shape wrt input
    else
      m1 =  (u*s1(-,))(,+) * vt(+,); // output has same shape as input
    return m1;
}

func check_eigen(l, &nmf, msg, disp=, thres=, condi=, win=,verb=)
// checking eigen values
{
  if( is_void(win) ) win=0;
  if( !is_void(condi) ) {
    thres = max(l)/condi;
    nmf = -1;
  }
  if( thres==[] ) thres=1e-4;
  if ( nmf==-1 ) {
    nn = where(l<thres);
    if(is_array(nn)) {
      nmf = numberof(l)-min(nn)+1;
    } else {
      nmf = 0;
    }
  }
  if( disp==1) {
    winorig = window();
    window, win;
    fma; limits; logxy,0,1;
    plg, l;
    plg,[l(1),l(0)],numberof(l)-[nmf,nmf],marks=0,color="red";
    if( msg==[] ) msg=" ";
    pltitle, msg;
    window, winorig;
  }
  if(verb)
    write,format="Number of filtered modes = %d \n",nmf;
}

func invert_eigen( L, nfilt )
// invert eigen values vector
{
    if( nfilt>0 ) {
        L1 = L;
        L1(1-nfilt:0) = 1.0;
        L1 = 1./L1;
        L1(1-nfilt:0) = 0.0;
    } else {
        L1 = 1.0 / L;
    }
    return L1;
}


/*
 _____ ____  ____   ___  ____  
| ____|  _ \|  _ \ / _ \|  _ \ 
|  _| | |_) | |_) | | | | |_) |
| |___|  _ <|  _ <| |_| |  _ < 
|_____|_| \_\_| \_\\___/|_| \_\
*/
func propagateError2Dm(eps, mca, lib=,verb= )
/* DOCUMENT cvv = propagateError2Dm(eps, mca, lib= )
   
     Computes the cov matrix of the residual error.
     Returns volts^2 when mca=mca is properly defined.
     
     The input covariance matrix <eps> is in arcsec^2.
     The output covariance matrix <Cvv> is in volts^2, to be used
     with functionfuncInflu(x,y,x0) that returns values in meters.

     SEE ALSO:
 */
{
  if(verb)
    write,format="Computing cov matrix on voltages (%s)\n",lib;
 
  Cvv = mca(,+) * (eps(,+)*mca(,+))(+,);
  
  // Cvv must be scaled to get volts^2 units, when the mca is in rd/volt.
  // When mca is in arcsec/volt, no unit change is required.
  if(0)
    Cvv /= 206265.^2;
  
  return Cvv;
}

func getMap( Caa, average )
/* DOCUMENT
     
   SEE ALSO:
 */

{
  extern data;
  
  // indices de la carte 2D carree de taille nsspXnssp, ou il y a 'vraiment' des sous-pupilles
  nn = *data.dm.csI;
  
  // creation du tableau des decalages
  np = data.dm.ndiam;
  xx = indgen(np)(,-:1:np)(nn); // coord en x des ssp
  dx = xx(,-) - xx(-,); // Ecart sur x entre chaque sous pupille
  yy = indgen(np)(-:1:np,)(nn); // coord en y des ssp
  dy = yy(,-) - yy(-,); // Ecart sur y entre chaque sous pupille

  // transformation des decalages en indice de tableau
  dx += np; 
  dy += np; 

  // transformation d'un couple de decalages (dx,dy) en un indice du tableau 'Map'
  Map = div = array(0.0, np*2-1, np*2-1);
  ind = dx(*)+(np*2-1)*(dy(*)-1);
  for(i=1;i<=numberof(ind);i++) {
    Map(ind(i)) += Caa(*)(i);
    div(ind(i)) += 1;
  }
  div( where(div==0) ) = 1;  // pour eviter division par 0
  if( average==1 )
    Map /= div;

  return Map;
}

func intersample( Cvvmap, N, ud, D, dactupix, couplingfact)
/* DOCUMENT res = intersample( Cvvmap, N, ud, D, dactupix, couplingfact)
   
   Cvvmap is the 'map' of the Cvv matrix (cov matrix of tomo error
   expressed on volts). The "volts" unit must be used together with
   the influence function funcInflu(x,y,dm.x0) expressed in meters.

   Then, the result of intersample is in meter^2.

   <ud>       : size of pixel of Dphi space (meters)
   <dactupix> : inter-actuator pitch expressed in pixels of size ud
   <couplingfact> : coupling factor of one influ func to the neighbour
   SEE ALSO:
 */
{
  if(verb)
    write,format="%s ... ","Interpolating Dphi map"

  // computation of influence function
  x = (indgen(N)-(N/2+1)) * ud / (D/2.);   // x exprime en rayon pupille
  x = x(,-:1:N);
  y = transpose(x);
  fi = funcInflu(x,y,couplingfact);  // expressed in meters
  r = abs(x,y);
  pup = r<=1.0 & r>data.tel.obs;
  
  map = array(0.0,N,N);
  ncmap = N/2+1;
  
  // size of the side of Cvvmap (always odd number)
  ncov = dimsof(Cvvmap)(0);
  // nber of elements on each side of the center of Cvvmap
  nelem = (ncov-1)/2;
  rr = ncmap-nelem*dactupix:ncmap+nelem*dactupix:dactupix;
  map(rr,rr) = Cvvmap;
  if(verb)
    write,"done";
  
  // Computing the phase correlation function
  // One should have corr(0) = phase_variance.
  // Computing corr(0) is done using the <v_i^2> (diagonal of Cvv).
  // We decided that <v^2> is the average value for 1 single actuator (i.e. it's the average
  // of the diagonal of Cvv).
  // Then the average phase variance over the pupil equals to
  // (1/S_pupil) * $_pupil(fi^2) * Nactu * <v^2>
  // with S_pupil the surface, and $ is an integral. Nactu needs to be here
  // because if it wasn't, we'd have computed the phase variance over the pupil
  // with only 1 single actu moving.
  // So, in our formula, we have replaced the value of (S_pupil/Nactu) by (dactu^2).
  // The (dactu^2) needs to be expressed in pixels because our integral $(fi^2) is not
  // a real integral : it's just summing pixels instead.
  corr = fft(abs(fft(fi))^2 * fft(map), -1).re / (numberof(fi) * dactupix^2) ;

  // From correlation to Dphi
  // Dphi(r) = 2*C(0) - 2*C(r)
  // We take advantage we need to do a multiplication to multiply by another factor
  // in the same line. This is to translate dphi from m^2 into rd^2
  fact = 2 * (2*pi/data.camir.lambda_ir)^2;
  dphi = fact*corr(ncmap,ncmap) - fact*corr;
  dphi = roll(dphi);

  // and ... we don't do that any more
  //  dphi1  = dphi / (abs(FTOtel)+1);
  
  return dphi;

  // computation of the PSF
  pup = FTOtel > FTOtel(1,1)/1e9;
  k = 2*pi/data.lambda.camir;
  psf = roll(fft(exp(-0.5 * k * k * dphi * pup) * FTOtel).re);

  
  return psf;
}


func dphiFromUij(Cvv,mode,verb=)
{
  tic,3;
  if(verb)
    write,format="%s ... ","Interpolating Dphi map";
  if(is_void(mode)) mode = "Vii";
  
  // computation of influence function and defining the pupil
  N = data.fourier.npix;
  x = (indgen(N)-(N/2+1)) * data.fourier.ud / (data.tel.diam/2.);   // x exprime en rayon pupille
  x = x(,-:1:N);
  y = transpose(x);
  P = circularPupFunction(x,y,data.tel.obs);

  //autocorrelation of the pupil
  fftpup = fft(P);
  conjpupfft = conjugate(fftpup);
  G = fft(fftpup*conjpupfft,-1).re +1e-15;
  
  
  //defining the modes
  if(mode != "intersample"){
    fi = readfits("dm_modes.fits",err=1);
    if(is_void(fi)){
      fi = array(0.,N,N,data.dm.nactu);
      xi = (*data.dm.csX);
      yi = (*data.dm.csY);
      for(i=1;i<=data.dm.nactu;i++){
        fi(,,i) = funcInflu(x-xi(i),y-yi(i),data.dm.x0);
      }
      writefits,"dm_modes.fits",fi*double(P(,,-));
    }
  }else{
    fi = funcInflu(x,y,data.dm.x0);
  }

  if(mode == "Uij"){
    dphi = 0*P;
    for(i=1;i<=data.dm.nactu;i++){
      for(j=i;j<=data.dm.nactu;j++){
        //modes
        Mi = fi(,,i);
        Mj = fi(,,j);
        //Uij computation
        t1 = (fft(Mi*Mj)*conjpupfft).re;
        t2 = (fft(Mi)*conjugate(fft(Mj))).re;
        Uij = fft(2*(t1 - t2),-1);
        //summing modes
        dphi +=  ((i!=j)+1.)*Cvv(i,j)*Uij;
      }
    }
      //multiplynig by a mask
      dphi = dphi.re/G;
      G /= max(G);
      msk = G > 1e-5;
      dphi *= msk/2.;
  }else if(mode == "Vii"){
    //Diagonalizing the Cvv matrix
    l = SVdec(Cvv,Bt,B);
    M = fi(,,+)*B(,+);
  
    //loop on actuators
    tmp = Mi = 0.*P;
    for(i=1;i<=data.dm.nactu;i++){
      Mi =  M(,,i);
      //Vii computation
      Vii = (fft(Mi*Mi)*conjpupfft).re - abs(fft(Mi))^2;
      //summing modes into dm basis
      tmp += l(i) * Vii;
    }
    dphi = fft(tmp,-1).re/G;
    //multiplynig by a mask
    msk = G/max(G) > 1e-5;
    dphi *= msk;
  
  }else if(mode == "intersample"){
    map = array(0.0,N,N);
    dactupix = data.fourier.dactupix;
    Cvvmap = getMap(Cvv,1);
    ncov = dimsof(Cvvmap)(0);
    // nber of elements on each side of the center of Cvvmap
    ncmap = N/2+1;
    nelem = (ncov-1)/2;
    rr = ncmap-nelem*dactupix:ncmap+nelem*dactupix:dactupix;
    map(rr,rr) = Cvvmap;
    corr = fft(abs(fft(fi))^2 * fft(map), -1).re / (numberof(fi) * dactupix^2) ;
    dphi = corr(ncmap,ncmap) - corr;
    dphi = roll(dphi);
  }

  
  //getting the phase structure function in rd^2
  dphi *= 2*(2*pi/data.camir.lambda_ir)^2.;

  if(verb)
    write,format="done in %.3g s\n",tac(3);

  return dphi;
}

/*
 ____ ___ __  __ _   _ 
/ ___|_ _|  \/  | | | |
\___ \| || |\/| | | | |
 ___) | || |  | | |_| |
|____/___|_|  |_|\___/ 
*/



func generateRandomVolts(cvv, nsample, verb)
{
  if(verb)
    write,format="Generating %d random draws on %d actuators ....... ",nsample,dimsof(cvv)(0);
  nactu = dimsof(cvv)(0);
  ll = SVdec(cvv, u);
  u *= sqrt(ll)(-,);

  ran = u(,+) * random_n(nactu, nsample)(+,);
  if(verb)
    write,"done";
  return ran;
}

func generatePsf( ran, N, champ_Dphi, D, couplingfact,verb= )
{
  extern dm;
  
  

  if(verb)
    write,format="%s ... ","Interpolating Dphi map"
  ud = champ_Dphi/double(N);  // taille des pixels dans espace Dphi (metres)

  pitch = D / (dm.Ndiam-1);   // pitch actuators in meters
  pitch /= ud;                // pitch actuators in pixels
  pitch = long(pitch+0.01);   // pitch actuators in pixels (round to integer value)

  // computation of influence function
  x = (indgen(N)-(N/2+1)) * ud / (D/2.);   // x exprime en rayon pupille
  x = x(,-:1:N);
  y = transpose(x);
  fi = funcInflu(x,y,couplingfact);  // expressed in meters
  r = abs(x,y);
  pup = r<=1.0 & r>tomo.tel.obs;
  
  map = array(0.0,N,N);
  ncmap = N/2+1;
  
  // size of the side of Cvvmap (always odd number)
  ncov = dimsof(Cvvmap)(0);
  // nber of elements on each side of the center of Cvvmap
  nelem = (ncov-1)/2;
  rr = ncmap-nelem*pitch:ncmap+nelem*pitch:pitch;
  map(rr,rr) = Cvvmap;
  if(verb)
    write,"done";
  
  FTOtel = fft(abs(fft(pup))^2, -1).re;
  dphi = fft(abs(fft(fi))^2 * fft(map), -1).re / roll(abs(FTOtel)+10000);
  dphi = dphi(ncmap,ncmap) - dphi;
  dphi = roll(dphi);
  
  // in fact dphi is not dphi, it should be multiplied by 2 to be a
  // real dphi, but one needs to multiply it by 1/2 in the
  // exp(-1/2.dphi)
  map = array(0.0,N,N);
  ncmap = N/2+1;
  
  // size of the side of Cvvmap (always odd number)
  ncov = dimsof(Cvvmap)(0);
  // nber of elements on each side of the center of Cvvmap
  nelem = (ncov-1)/2;
  rr = ncmap-nelem*pitch:ncmap+nelem*pitch:pitch;
  map(rr,rr) = Cvvmap;
  if(verb)
    write,"done";
  
  FTOtel = fft(abs(fft(pup))^2, -1).re;
  dphi = fft(abs(fft(fi))^2 * fft(map), -1).re / roll(abs(FTOtel)+10000);
  dphi = dphi(ncmap,ncmap) - dphi;
  dphi = roll(dphi);
  
  // in fact dphi is not dphi, it is only 0.5*dphi but we keep it like
  // this as we need to multiply it by 1/2 anyway in the
  // exp(-1/2.dphi)
  //
  pup = FTOtel > FTOtel(1,1)/1e9;
  dphi *= pup;
  k = -(2*pi/1.65e-6)^2.;
  psf = roll(fft(exp(k * dphi) * FTOtel).re);
  
  return psf;
}









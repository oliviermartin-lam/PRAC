func computeOTFtomographic(averageMode,&dphi,para=,Cee=,verb=)
/*DOCUMENT dphiFromUij(Cvv,"Vii")

  Returns the phase structure function (SF) average on the pupil from the tomographic error in rd^2. You have to specified the inputs mode as either "Uij" to average the  SF using Veran's method (1997), or "Vii" for using Gendron's method (2006) or "intersample" to use simplified Gendron's approach (2014 Kaust).
  The covariance matrix of voltages Cvv must be given in volts^2.

 */
{
  tic,3;

  if(is_void(Cee))
    Cee = computesCeeMatrix(rtc.obsMode,para=para,verb=verb);

  //computes the covariance matrix of voltages in volts^2
  Cvv = propagateError2Dm(Cee, *rtc.mc );
  
  if(verb)
    write,format="%s ... ","Interpolating Dphi map";
  if(is_void(averageMode)) averageMode = "Vii";
  
  // computation of influence function and defining the pupil
  N = tel.nPixels;
  P = circularPupFunction(tel.diam,tel.obs,tel.nPixels,tel.pixSize) ;
  //pupil function expressed in pupil radius

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
      OTF_tomo  = exp(-0.5*dphi);

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
    OTF_tomo  = exp(-0.5*dphi);

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
    OTF_tomo  = exp(-0.5*dphi);
  }

  if(verb)
    write,format="done in %.3g s\n",tac(3);

  
  return OTF_tomo;
}



    
  
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
  its = rtc.its;
  
  // Creation of subaperture points of the TS
  nssp = wfs(its).nLenslet;
  xy =  compute_ri(nssp, 2.0, tel.obs,0,0);

  // matrix of distances (in x and y) from subaperture centre to DM
  dx = xy(,1) - (*dm.csX)(-,);
  dy = xy(,2) - (*dm.csY)(-,);

  // subap size
  ssize = 2./nssp;
  
  // creation of MI
  nsubap = int(wfs(its).nValidSubAp);
  mia = array(0.0, nsubap*2, dm.nActu);

  // here <mia> has no unit, it is just a difference between the edges of a subap
  // The result is expressed in meters.
  mia(1:nsubap,)  = difffuncInflu(dx, dy, ssize, dm.x0);
  mia(nsubap+1:,) = difffuncInflu(dy, dx, ssize, dm.x0);

  // here we convert mia into an angle (in rad) per volt
  // We need to divide the OPD by the size of subap to get an angle (radians)
  mia /= tel.diam/nssp;
  //we switch from radians to arcseconds
  mia *= radian2arcsec;
  
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
  
func invgen(mat , nfilt,cond=, t=)
/* DOCUMENT invgen( mat, nfilt, t=1 )
	mat   : matrix to be inverted
	nfilt : number of modes to be filtered
        t=1   : output has transposed shape wrt input
        t=0 (default) : output has same shape as input
*/
{
    s = SVdec(mat,u,vt);
    if(cond){
      thres = max(s)/cond;
      nn = where(s<thres);
      if(is_array(nn))
        nfilt = numberof(s)-min(nn)+1;
      else
        nfilt = 0;
    }
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
func propagateError2Dm(eps, mca)
/* DOCUMENT cvv = propagateError2Dm(eps, mca, lib= )
   
     Computes the cov matrix of the residual error.
     Returns volts^2 when mca is properly defined.
     
     The input covariance matrix <eps> is in arcsec^2.
     The output covariance matrix <Cvv> is in volts^2, to be used
     with functionfuncInflu(x,y,x0) that returns values in meters.

     SEE ALSO:
 */
{
 
  Cvv = mca(,+) * (eps(,+)*mca(,+))(+,);
  
  // Cvv must be scaled to get volts^2 units, when the mca is in rd/volt.
  // When mca is in arcsec/volt, no unit change is required.
  if(0)
    Cvv /= radian2arcsec^2;
  
  return Cvv;
}

func getMap( Caa, average )
/* DOCUMENT
     
   SEE ALSO:
 */

{
  // indices de la carte 2D carree de taille nsspXnssp, ou il y a 'vraiment' des sous-pupilles
  nn = *dm.csI;
  
  // creation du tableau des decalages
  np = dm.nActu;
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





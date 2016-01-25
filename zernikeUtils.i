
func nmOrder(i, &n, &m)
/* DOCUMENT nm = ordreNM(i, n, m)
     Pour un zernike donne, defini par son indice i, calcule
     l'ordre azimutal m et l'ordre radial n.
 */
{
  n = int( (-1.+sqrt(8*(i-1)+1))/2.);
  p = (i-(n*(n+1))/2);
  k = n%2;
  m = int((p+k)/2)*2 - k;
}


func varZernike( s , icam ,arc=)
/* DOCUMENT varzer = varZernike( s , icam )
     Compute variance of zernike polynomials.
     The variable <s> can either be a slope vector, or an array of slopes.
     For a single vector, the square of the zernike coeffs is returned.
     For a series, the centered mean-square value is returned.

     The variance array returned goes from TILT to Z35
   SEE ALSO:
 */
{ 
  // zernike reconstruction on zernike modes
  z = (*sys.slopesToZernikeMatrix)(,+)*s(+,);
  // variance
  if(dimsof(s)(1)==1 || dimsof(s)(3)==1) {        // if s has 1 slope
    varz = z^2.;
  } else {
    varz = z(,rms)^2.;
  }
  units = wfs(icam).unit*1000.;
  if(arc) units = 0.5*tel.diam*1e9/radian2arcsec;
  varz *= units^2;      // in nm^2
  return varz;
}

func zernikeSpectrum(radorder,dr0,dl0)
/* DOCUMENT zernikeVonKarmanSpectrum(jmax, dlo, dro)
   
     <jmax> is the number of the last Zernike,
     <dlo> is D/L0
     <dro> is D/r0
     The function returns the theoretical spectrum of Zernike coeffs,
     with a given r0 and outer scale. It uses a formula given by Rodolph Conan in his phD thesis.

     The returned vector ranges from Z_2 (tilt) to Z_jmax.
     The piston variance is not returned (it could, but it's not ...).

     The function works for large L0>D/4, i.e. for 0<= (D/L0) < 4,
     it will diverge for L0<D/4.

     Number of terms of the taylor expansion.
     This number can be increased if higher accuracy is required,
     e.g. for computing spectra for small L0<D/4, however the success
     is not guaranteed.
     
   SEE ALSO:
 */
{
  local n,m;

  jmax = ((radorder+1)*(radorder+2))/2;
  varZernike=array(0.0, jmax-1);

  zerfactor = [1.00127,1.00127,1.0089,1.00421,1.01433,1.05697,1.05697,1.0464,1.0464,1.07974,1.12548,1.03053,1.07697,1.05695,1.0647,1.0647,1.13497,1.13497,1.12401,1.12401, 1.46191,1.32155,1.30285,0.752463,1.25776,0.911428,0.994586,2.69474,2.69474, 2.18867,2.18867,1.3627,1.3627,1.39158,1.39158];

  factNorm = 2*gamma_R(11./6.)*((24./5)*gamma_R(6./5))^(5./6)/pi^(3./2.)*abs(dr0)^(5./3);

  nmOrder,indgen(2:jmax),n,m;
  nmax = max(n);

  pidf0 = pi*abs(dl0)
  if(pidf0>pi*4) {
    pidf0=pi*4;
  }
  
  kmax=50;
  for(i=1; i<=nmax; i++) {
    s = series_sum_diag(i,kmax,pidf0) * (i+1) * factNorm;
    varZernike(where(n==i)) = s;
  }

  varZernike *= zerfactor;

  return varZernike;
}



func series_sum_diag(n,kmax,pidf0)
/* DOCUMENT series_sum_diag(n,kmax,pidf0)

   This function computes a sum of terms composed with Gamma
   functions, for the computation of zernike variances with finite
   L0. The formula is extracted from R. Conan thesis (equation 4.16
   page 116).
   The formula has been modified to use the relation
   Gamma(k+x)=(k-1+x).Gamma(k-1+x) and thus avoid to compute the Gamma
   functions at each iteration.
   
   SEE ALSO:
 */
{
  n2 = 2*n;

  // compute n! = Gamma(1+n)
  fn = 1.00;
  // compute (2+n1+n2)! using the previous results again
  for(i=2; i<=2+n2; i++)
    fn*=i;  // (2+n+n)! = Gamma(3+n+n)

  // computation of all gamma_R as a whole
  pregamma = gamma_R([1.5 + n, 5./6 - n, n-5./6, n+23./6]);

  // Initialisation
  uk = pregamma(1) * pregamma(2) / ((n+1) * fn);  // u0

  vk = pregamma(3)*gamma_R(7/3.)*gamma_R(11./6)/gamma_R(17./6) / pregamma(4);
    
  pidf0_n253 = pidf0^(n2-5/3.);
  pidf0_2 = pidf0*pidf0;
  pidf0_2k = 1.00;   // preparing pidf0^(2*k), for k=0
  fk = 1.0;          // preparing k!, for k=0

  s = (uk * pidf0_n253 + vk);    // s(k=0)
  
  for(k=1; k<=kmax; k++) {
    fk *= k;    // k!
    pidf0_2k *= pidf0_2;   // computation of  pidf0^(2*k)
    
    uk *= ((0.5+n+k)*(k+n))/((2+n2+k)*(1+n+k)*(5./6-n-k));
    vk *= (k+4./3)*(k+5./6)/((n-5./6-k)*(n+17./6+k)*(11./6+k));
                                       
    tmps = (pidf0_n253 * uk + vk) * pidf0_2k / fk;
      
    if(k%2) s -= tmps;
    else  s += tmps;                                     
  }
  
  return s;
}

func computeZernikePolynomials(r,t,i)
  /* DOCUMENT RTzer(r,t,i)

  Meme chose que zer(r,t,i), mais ici r et t sont des tableaux 1D,
  et le produit externe est effectue dans la fonction.
  Du coup, le polynome(r) n'est calcule que sur les 1D-points de r,
  le sin ou cos(m.theta) n'est calcule que sur les 1D-points de theta,
  seul le produit est fait en externe.

  Du coup, ca va TRES vite.

  Pour afficher :
  xx = r(,-)*sin(th)(-,);
  yy = r(,-)*cos(th)(-,);
  plf, RTzer(r,t,52),xx,yy;
  
   */
{
  if(i==1) {
    return 0*r+1.;
  }
  // calcul de n et m a partir de i
  nmOrder,i, n, m
  a = polyfute(m,n);
  
  Zi = computePolynome(n,m,a,r) * sqrt(n+1) * sqrt(2);
  
  if( m!=0 ) {
    if( i%2 ) {
      Zi = Zi * sin(m*t);
    }
    else {
      Zi = Zi * cos(m*t);
    }
  }
  
  return Zi;
}

func computePolynome(n,m,a,r)
  /* DOCUMENT evaluate_poly(n,m,a,r)
     n is the radial order
     m is the azimutal order
     a[] is the list of coefficient, with a(i+1) the coeff of r^i
     r is the variable of the polynomial
  */
{
  if( n>1 )
    r2 = r*r;

  p = a(n+1);
  for(i=n-2;i>=m;i-=2) {
    p = p*r2 + a(i+1);
  }
  
  if(m==0) return p;
  else if(m==1) p*=r;
  else if(m==2) p*=r2;
  else p*=r^m;
  
  return p;
}




func polyfute(m,n)
  /* DOCUMENT polyfute(m,n)

  Les coefs des poly de zer sont des K_mn(s).
  Le coeff K_mn(s) pondère r^(n-2s)

  Il y a la relation de recurrence
  K_mn(s+1) =  K_mn(s) * ((n+m)/2-s)*((n-m)/2-s)/(s+1)/(n-s)

  Il y a aussi
  K_mn(0) = n! / ((n+m)/2)! / ((n-m)/2)!
  
  */
{
  a = array(double,n+1);

  // Calcul de K_mn(0)
  st = 2;                           // start index for dividing by ((n-m)/2)!
  coef = 1.00;
  for(i=(n+m)/2+1; i<=n; i++) {     // calcul de  n! / ((n+m)/2)!
    if( st<=((n-m)/2) & i%st==0 ) {
      j = i/st;
      st++;
      coef *= j;
    } else {
      coef *= i;
    }
  }
  // division by ((n-m)/2)! (has already been partially done)
  for(i=st;i<=(n-m)/2;i++) coef /= i;

  a(n+1) = floor( coef + 0.5 );   // pour K_nm(0)
  
  for(i=1;i<=(n-m)/2;i++) {
    coef *= -((n+m)/2-i+1)*((n-m)/2-i+1);
    coef /= i;
    coef /= (n-i+1);
    a(n-2*i+1) = floor( coef + 0.5 );
  }
  return a;
}

func fitKolmo(x,a)
{
  radorder = x(1);
  i0 = x(2);
  i1 = x(3);
  tmp = zernikeSpectrum(radorder, a(1), a(2))(i0:i1);

  return log(abs(tmp));
  
  // abs is not required, but allows to avoid bugs when running on crap data ..

}

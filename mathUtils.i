func arrondi(a,n)
/* DOCUMENT b = arrondi(a,n)
   Retourne l'arrondi le nombre n-decimal le plus proche de a.
*/
{

  
  up = ceil(a*10.^n)/10.^n;
  down =  int(a*10.^n)/10.^n;

  if(numberof(a)==1){
    
    if(abs(a - up) <= abs(a-down)) return up;
    else  return down;
  }
  
  else{
    na = numberof(a);
    res = array(0.,na);
    
    for(k=1;k<=na;k++){
      if(abs(a(k) - up(k)) <= abs(a(k)-down(k))){ res(k) = up(k);}
      else{ res(k) = down(k);}
     
    }
    return res;
  }
  
}
func autocorrelation( a )
/* DOCUMENT aa = autocorrelation( a )

   computes the autocorrelation so that
   max(aa) == sum(a^2)
     
   SEE ALSO:
 */
{
  b = abs(fft(a));
  b = fft(b*b).re;
  n2 = numberof(a);  // N*N
  b /= n2;
  return b;
}
func determinesAverageOnProfile(cn2h,q,qmin,qmax)
{

  if(!anyof(q))
    return 0;
  q2 = q;cn2h2 = abs(cn2h);
  //truncation of the outliers layers
  
    if(is_void(qmin))
      qmin = min(q2);
    if(is_void(qmax))
      qmax = max(q2);
    w = where(q2>= qmin & q2<= qmax);

    if(is_array(w)){
      q2 = q2(w);
      cn2h2 = cn2h2(w);
    }
  
  tmp = (sign(q2)*cn2h2*abs(q2)^(5/3.))(sum);
  qavg = sign(tmp)*(abs(tmp)/(cn2h2(sum)))^(3/5.);
 
  return qavg;

}
func computesIsoplaneticAngleFromProfile(cn2h,altitude,rzero=)
{
  havg = determinesAverageOnProfile(cn2h,altitude);
  if(is_void(rzero))
    rzero = sum(cn2h)^(-3/5.);

  if(havg !=0)
    theta0 = (0.314*rzero/havg)*206265.;
  else
    theta0 = 0;
  return theta0;
}
func computeModesToBeFiltered(Caa,TScam,listWfsType=,condCaa=)
/* DOCUMENT   nmodes = computeModesToBeFiltered(Caa,listWfsType=);

   Determines the optimal number to filter directly from the Caa
   
   Olivier Martin.
 */
{

  if(is_void(listWfsType)) listWfsType = wfs.type;
  if(is_void(TScam)) TScam = rtc.its;
  // SVD decomposition
  eigenval = SVdec(Caa);
  // .....number of modes to be kept.....//

  // 72 for WFS of type 1, 70 for WFS of type 2 (laser), 2 for WFS of type==3 (tilt only)
  modes_garder1 = sum((wfs(others(rtc.its)).type==1)*(wfs(others(rtc.its)).nValidMeas)) +
    sum((wfs(others(rtc.its)).type==2)*(wfs(others(rtc.its)).nValidMeas-2)) +
    sum((wfs(others(rtc.its)).type==3)*(2));
  
  eigenval_type = eigenval(1:modes_garder1);
  if( min(eigenval_type)==0 )
    eigenval_type += 1e-100;
  
  // spots where the eigenvalue curve "jumps"
  // the jumps that occur in the 10 first values are not considered as
  // a real "drop"
  seuil=1;
  tmp = (-log(eigenval_type)(dif)>seuil);
  ntmp = 10;

  if(numberof(tmp)>ntmp){jumps = where(tmp(ntmp:));}
  
  else{
    jumps = where(tmp);
    ntmp=1;
  }

  if( !is_array(jumps) ) {   // if no 'jump'
    modes_garder2 = modes_garder1;
  }
  else {
    modes_garder2 = jumps(1)+ntmp-1; // the first 'jump' (after the first 10) determines where we cut
  }

  modes_garder = min(modes_garder1, modes_garder2);

  if(condCaa){
    lmin = max(eigenval_type)/condCaa;
    w = where(eigenval_type<lmin);
    if(is_array(w))
      modes_garder = w(1)-1;
  }
  //number of modes to be filtered = total modes - TSmodes - typemodes
  nmodes = numberof(eigenval(norange(rtc.its))) - modes_garder;
    
  return nmodes;
}
func conjugate(c)
/* DOCUMENT cstar  =conj(c)
   Returns the conjugate value of the complex number c;
 */
{
  if(is_complex(c)){
    return c.re - c.im;
  }else{
    return c;
  }
}
func createPixarcVector(void)
/* DOCUMENT pixarc = createPixarcVector()

   Creates a vector of length data.Nslopes containing the pixarc of each slope.
   
   SEE ALSO:
 */
{
  pixarc = array(1., rtc.nSlopes);
  for(i=1; i<=rtc.nWfs; i++) {
    pixarc(slrange(i)) = wfs(i).pixSize;
  }
  return pixarc;
}
func identite(n)
{
m = array(0.0,n,n)
m(*)(1::n+1) = 1
return m
}
func inverse_mat(mat , nfilt, &s)
/* DOCUMENT invgen( mat, nfilt, s)
	<mat>   : matrice a inverser
	<nfilt> : nombre de modes a filtrer
        <s>     : output parameter: eigenvalues
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
    m1 =  (u*s1(-,))(,+) * vt(+,);
    return m1;
}
func getFWHM(g,Fe)
{

  nf = numberof(g);
  // determining the FWHM on the y autocorrelation
  npt =nf/2;
  k=2;
  while( g(k)>0.5 & k<npt) k++;   // search first point below 0.5
  // linear regression on 5 points around the half-max point
  if( (k-2)>3 & (k+2)<npt ) {
    y = g(k-2:k+2);
    x = indgen(k-2:k+2);
    tmp = (avg(x*x)-avg(x)^2.);
    a = (avg(y*x)-avg(x)*avg(y))/tmp;
    b = avg(y) - a*avg(x);
    k = (0.5-b)/a;  // solve equation a*k+b = 0.5
  }
  tau0 = k /Fe;   // tau0 in seconds

  return tau0;
}


func giveTomoMode(wfstype,OBS_MODE,RECTYPE)
/*DOCUMENT
 */
{

  NLGS = sum(wfstype==2);
  NNGS = sum(wfstype==1) - wfstype(0);//subtract TS
  NTTGS = sum(wfstype==3);

  isLGS = isNGS = isTTGS ="";
  if(NLGS) isLGS = var2str(NLGS) + "L";
  if(NNGS) isNGS = var2str(NNGS) + "N";
  if(NTTGS) isTTGS = var2str(NTTGS) +"T";

  if(RECTYPE == "glao" | RECTYPE == "GLAO") OBS_MODE = "GLAO";
  if(OBS_MODE == "SCAO") return OBS_MODE;
  
  TomoMode = OBS_MODE + isLGS + isNGS + isTTGS;

  return TomoMode;
}
func minmax(arg) 
/* DOCUMENT minmax(arg)
 * Returns a vector containing the min and the max of the argument
 * F.Rigaut 2001/09
 * SEE ALSO:
 */
{
  return [min(arg),max(arg)];
}
func proportionalRegress(y,x)
/*DOCUMENT alpha = proportionalRegress(y,x)
  
  Returns the proportionnal regression coefficient alpha which
  minimizes sum((y - a*x)^2 ).
  
*/

{
  x = double(x);
  y = double(y);
  return sum(y*x)/sum(x^2);
}
func regress(y,x,pen=,ab=)
/* DOCUMENT regress(y,x)
            regress,y,x,pen=1
            coef = regress(y,x,ab=1)
   returns the regression coefficient,
   returns [a,b] coeffs when ab=1
   Prints A and B (y=Ax+B) when pen=1
 */
{
  tmp = (avg(x*x)-avg(x)^2.);
  a = 0.0;
  if( tmp!=0 )
    a = (avg(y*x)-avg(x)*avg(y))/tmp;
  b = avg(y) - a*avg(x);
  if( !is_void(pen) ) {
    print,"pente = ",a;
    print,"orig  = ",b;
  }
  if( !is_void(ab) ) {
    return [a,b];
  }
  return a;
}
func SQRT(x)
{
  if(x==[]) return 0.;
  return sqrt(abs(x)) * sign(x);
}
func unbiasesMtFromSensitivities(mt)
/*DOCUMENT

 */
{
  mt2 = mt;
  pixarc =  createPixarcVector();
  mt2/= pixarc(norange(rtc.its))(-,);
  mt2*= pixarc(slrange(rtc.its));

  return mt2;
}
func tic(counterNumber)
/* DOCUMENT tic(counter_number)
 * Marks the beginning of a time lapse
 * ex: tic ; do_something ; tac()
 * will print out the time ellapsed between tic and tac
 * a counter number can optionaly be specified if several
 * counters have to be used in parallel.
 * F.Rigaut 2001/10
 * SEE ALSO: tac
 */
{
  if (counterNumber == []) counterNumber = 1;
  if (counterNumber > 10) error,"tic and tac are limited to 10 time counters !";

  el = array(double,3);
  timer,el;
  _nowtime(counterNumber) = el(3);
}
if (numberof(_nowtime)!=10) _nowtime = array(double,10);
func tac(counterNumber)
/* DOCUMENT tac(counter_number)
 * Marks the end of a time lapse
 * ex: tic ; do_something ; tac()
 * will print out the time ellapsed between tic and tac
 * a counter number can optionaly be specified if several
 * counters have to be used in parallel.
 * F.Rigaut 2001/10
 * SEE ALSO: tic
 */
{
  if (counterNumber == []) counterNumber = 1;

  el = array(double,3);
  timer,el;
  elapsed = el(3)-_nowtime(counterNumber);

  return elapsed;
}

func takesDiag(Mat)
/* DOCUMENT diag(matrix)

Returns the diagonal elements of a matrix.
Fails for rectangular matrices.
*/
{
  n = dimsof(Mat)(2);
  return Mat(*)(1::n+1);
}
func where4Vector(w1,w2,notequal=)
/*DOCUMENT
 */
{
    
  n1 = numberof(w1);
  n2 = numberof(w2);

  wmin = w1;wmax = w2;
  if(n1>=n2){
    wmin = w2;wmax = w1;
  }
  n = numberof(wmax);
  wVect =  wVectnot = array(0.,n);
   
  for(i=1;i<=n;i++){
    tmp = sum((wmax(i) == wmin));
    if(tmp)
      wVect(i) = tmp;
    if(notequal && !tmp)
      wVectnot(i) = 1;
  }
  if(notequal)
    return wVectnot;
  else
    return wVect;
}






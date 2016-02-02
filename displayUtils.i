/*
 ____  _           _                __              _                          
|  _ \(_)___ _ __ | | __ _ _   _   / _| ___  _ __  | |    ___  __ _ _ __ _ __  
| | | | / __| '_ \| |/ _` | | | | | |_ / _ \| '__| | |   / _ \/ _` | '__| '_ \ 
| |_| | \__ \ |_) | | (_| | |_| | |  _| (_) | |    | |__|  __/ (_| | |  | | | |
|____/|_|___/ .__/|_|\__,_|\__, | |_|  \___/|_|    |_____\___|\__,_|_|  |_| |_|
            |_|            |___/                                               
*/
func displayLayers(cn2h,alt,xylabels,l0=,percent=,col=, thick=,strength_rms=)
  /* DOCUMENT displayLayers,data.learn.cnh,data.learn.altitude

     Displays the turbulent profile with strengths and altitues given
     by cn2h and alt. If percent=1, plots turbulent profile with
     relative strength in percent. If strength_rms is an array of the
     same size than percent, it will plot black error bars on the
     strength.

 */
  
{

  // getting profile data from tomo structure
  x = cn2h;
  // if(!l0){
    //w = where(cn2h!=0);
    //x(w) =0.103/cn2h(w)^(-3/5.);
  //}
  //x = cn2h;
  if(!sum(sign(cn2h)==-1) == numberof(cn2h))
    x = abs(cn2h);
  if(percent){
    x=100*x/sum(abs(x));
    if(is_array(strength_rms))
      strength_rms = x*(strength_rms/abs(cn2h) - sqrt(sum(strength_rms^2))/sum(abs(cn2h)));
  }

  y = alt/1000.;
  nblayers = numberof(x);
  
  // organizing data
  nn = sort(y);
  y=y(nn);
  x=x(nn);

  // setting limits of display, so that the window is slightly larger than the plot, and
  // properly scaled when there s only 1 point to display
  ux = grow(x,0);
  uy = grow(y,0);
  dx = abs(ux(ptp)) * 1.1;
  dy = abs(uy(ptp)) * 1.1;

  //if( dx<1e-4 ) dx=1.0;
  if( dy<10 ) dy=1000.0;
  cx = (max(ux)+min(ux))/2.;
  cy = (max(uy)+min(uy))/2.;
  limits, cx-dx/2, cx+dx/2, cy-dy/2, cy+dy/2;
  // display blue bars .............
  for(i=1; i<=nblayers; i++) {
    if(is_void(col))
      col = [char(249)];   // dark blue
    if(is_void(thick))
      thick = dy/50.;

    plfp, col, [1,1,1,1]*y(i)+[-1,-1,1,1]*thick,[0,1,1,0]*x(i), [4],edges=0;
    
    if(!is_void(strength_rms)){
      plg, [y(i),y(i)], [x(i)-strength_rms(i)/2.,x(i)+strength_rms(i)/2.], marks=0, width=4;
      plg, y(i) + 1.1*thick*[-1,1],[x(i)-strength_rms(i)/2.,x(i)-strength_rms(i)/2.],marks=0,width=4;
      plg, y(i) + 1.1*thick*[-1,1],[x(i)+strength_rms(i)/2.,x(i)+strength_rms(i)/2.],marks=0,width=4;
    }
  }
  // display axes
  //plg, [0,0], [0,cx+dx/2], marks=0, width=5;
  //plg, [0,cy+dy/2], [0,0], marks=0, type=3;
  if(is_void(strength_rms)) strength_rms = 0*x;
  limits, cx-dx/2, cx + dx/2. + max(strength_rms/2.), cy-dy/2, cy+dy/2;

  if(is_void(xylabels)){
    if(percent){
      xytitles,"Relative contribution (%)","Altitude (km)";
    }else{
      xytitles,"Seeing at 500 nm  (arcsec) ","Altitude (km)";
    }
  }else{
    xytitles,xylabels(1),xylabels(2);

  }

}

func dispCoeffs(coeffs)
/* DOCUMENT Displays in the appropriate fields of the GUI
            the results of the learn process.
   
   SEE ALSO:
 */
{
  // L0, r0, altitude, noise, tracking, xaso, sensib and yaso are output variables retrieved from 'coeffs'
  unpackcoeffs, coeffs,  cnh, altitude, l0h;

  window,6;clr;limits,square=0;logxy,0,0;
  displayLayers, abs(cnh), altitude;
  window,7;clr;limits,square=0;logxy,0,0;

  w = where(l0h != 0);
  outscale = l0h;
  outscale(w) = abs(l0h(w));
  displayLayers, outscale, altitude,["Outer scale (m)","Altitude (km)"],l0=1;
  
}


func tracePupilles(coeffs, posx, posy)
{  
  // get parameters L0, r0, altitude
  unpackcoeffs, coeffs, cnh, alt,l0h;

  
  window, 8;clr;logxy,0,0;

  plt,"initial",0.58,0.85,color="red";
  plt,"fitted",0.58,0.83,color="green";

  // layer altitude on which pupils will be projected
  hmax = max(max(alt),100.);

  // drawing of the LGS printhrough on highest layer
  lgsH = wfs.lgsH;
  nn = where(lgsH>0);
  if( is_array(nn) )
    lgsH(nn) = 1./lgsH(nn);           // inverse of lgs altitude (0 for ngs)

  pupdiam = radian2arcsec*(tel.diam*max((1./hmax - lgsH),0));
  // rescaling of pupil size in terms of layer altitude and lgs altitude
  xaso = wfs.x;
  yaso = wfs.y;
  // fitted asterism
  fact = regress([xaso,yaso](*), [wfs.x,wfs.y](*),ab=1);
  fact = fact(1);
  fact = 1.00; // finalement, non. Enfin pour l instant.
  for(i=1;i<=rtc.nWfs;i++) {
    linetype = 1;                            // solid line
    if( i==rtc.its ) linetype = 3;          // dash-dot
    if( wfs.type(i)==0 ) linetype = 0;  // none
    cercle, pupdiam(i)/2, center=[xaso(i),yaso(i)], type=linetype, col=["blue","green","magenta","white"](wfs.type(i));
  }

  if(is_array(*rtc.ptrListLgs))
    plmk, yaso(*rtc.ptrListLgs)/fact, xaso(*rtc.ptrListLgs)/fact, color = "green", marker = 4, width=10, msize=0.8;

  if(is_array(*rtc.ptrListOffAxisNgs))
    plmk, yaso(*rtc.ptrListOffAxisNgs)/fact, xaso(*rtc.ptrListOffAxisNgs)/fact, color = "green", marker = 4, msize=0.8;
 
  
  // initial values that were initially in the tomo structure
  ix = wfs.x;
  iy = wfs.y;
 
  // theoretical asterism
  if(is_array(*rtc.ptrListLgs))
    plmk,iy(*rtc.ptrListLgs), ix(*rtc.ptrListLgs), color = "red", marker = 6, width=10;
  if(rtc.ptrListOffAxisNgs)
    plmk,iy(*rtc.ptrListOffAxisNgs), ix(*rtc.ptrListOffAxisNgs), color = "red", marker = 6;
  // central position
  plmk, posy, posx, color="blue", marker = 6;
  
  limits, square=1;
  limits;
}

func cercle(r, center=, col=,width=, type=)
/* DOCUMENT cercle([10,1,2])
   
     
   SEE ALSO:
 */
{
  if(is_void(center)){center = [0,0];}
  if(is_void(col)){col="red";}
  if(is_void(width)){width=1;}
  if(is_void(type)){type=1;}
  t = span(0,2*pi,111);
  plg, r*cos(t)+center(2), r*sin(t)+center(1), color=col, width=width, type=type ,marks=0;
}


/*
 _   _ _     _                                  
| | | (_)___| |_ ___   __ _ _ __ __ _ _ __ ___  
| |_| | / __| __/ _ \ / _` | '__/ _` | '_ ` _ \ 
|  _  | \__ \ || (_) | (_| | | | (_| | | | | | |
|_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_|
                      |___/                     
*/

func plotsBarDiagram(y,labs,y2=,y3=,thick=,col1=,col2=,col3=,title=,step=)
/*DOCUMENT
 */
{

  local y2,y3,thick,col1,col2,col3;
  if(is_void(step))
    step = 1;
  //number of bar to draw
  nbar = numberof(y);
  m = n = array(0.,nbar);
  fma;limits;
  if(is_void(thick))
    thick = 0.1;

  if(is_void(col1))
    col1 = [char(245)];   // dark blue
  if(is_void(col2))
    col2 = [char(251)]; //red
  if(is_void(col3))
    col3 = [char(250)]; //green
  
  
  if(!is_array(y2) && !is_array(y3) ){
    for(i=1; i<=nbar; i++){
      d=0;
      plfp, col1,[0,1,1,0]*y(i),([1,1,1,1]*i*step-d/2) + [-1,-1,1,1]*thick,[4];
      m(i) = max(y(i));n(i) = min(y(i));
    }
  }
  else if(!is_array(y3) ){
    for(j=1; j<=nbar; j++){
      if(step ==1)
        i = j;
      else
        i = j-1;
      d = 0.5*(0.5-thick);
      d = 2*thick;
      m(i) = max(y(j),y2(j));
      n(i) = min(y(j),y2(j));
      plfp, col1,[0,1,1,0]*y(j),([1,1,1,1]*i*step-d/2.) + [-1,-1,1,1]*thick,[4];
      plfp, col2,[0,1,1,0]*y2(j),([1,1,1,1]*i*step+d/2.) + [-1,-1,1,1]*thick,[4];
    }
  }
  else{
    for(i=1; i<=nbar; i++){
      d = 0.5*(0.5-thick);
      m(i) = max(y(i),y2(i));
      n(i) = min(y(i),y2(i));
      plfp, col1,[0,1,1,0]*y(i),([1,1,1,1]*i*step-d) + [-1,-1,1,1]*thick,[4];
      plfp, col2,[0,1,1,0]*y2(i),([1,1,1,1]*i*step) + [-1,-1,1,1]*thick,[4];
      plfp, col3,[0,1,1,0]*y3(i),([1,1,1,1]*i*step+d) + [-1,-1,1,1]*thick,[4];
      m(i) = max(y(i),y2(i),y3(i));
      n(i) = min(y(i),y2(i),y3(i));
    }
  }
  for(j=1;j<=nbar;j++){
    plt,labs(j),j - thick-d,max(abs(m(j))*1.05,10),tosys=1;
  }
  
  limits,0,step*nbar+1;
  if(min(n)>0)
    range,0,1.1*max(m);
  else
    range,1.1*min(n),1.1*max(m);
  if(title)
    xytitles," ","Wavefront error (nm rms)",[-0.01,0];  
}

func histo(t,n, color=, marks=,width=,type=)
/* DOCUMENT histo(tab,n)
        displays the histogram of array 'tab', on 'n' points.
        See also hist, or histogramme.
 */
{
  t = double(t);
  dt = abs(t(ptp));
  tmin = min(t) - dt/double(n);
  y = histogram( int(1+0.5+n*(t-tmin)/dt) );
  x = span(tmin,max(t),dimsof(y)(2));
  gplh, y,x,color=color,marks=marks,width=width,type=type;
}

func ihist(t, hcen, ihwidth, nice=)
/* DOCUMENT ihist(tab, hcen, ihwidth, nice=)
        displays the histogram of array 'tab', centered on 'hcen',
        and from/to hcen+/-ihwidth.
        The step of the histogram is 1.00.
        The parameter nice=1 allows to draw the avg and median values.
        
        See also hist, histo, or histogramme.
 */
{
  limits,square=0;
  fma; logxy,0,0;
  t = double(t);
  if( nice==1 ) {
    aaa = avg(t);
    mmm = median(t(*));
  }
  nn = where(abs(t-hcen)<=ihwidth);
  minimum = min(t(nn));
  y = histogram( t(nn) - minimum + 1 );
  x = indgen(dimsof(y)(2))+minimum-1;
  plg, y,x;
  if( nice==1 ) {
    plg,[0,max(y)],[mmm,mmm],marks=0,color="red";
    plg,[0,max(y)],[aaa,aaa],marks=0,color="blue";
  }
}

func histogramme(hmin,hmax,t,over=)
/* DOCUMENT histogramme(hmin,hmax,tab)
        displays the histogram of array 'tab', between 'hmin' and 'hmax'.
        See also histo, or histogramme.
 */
{
  hist(t(where((t>hmin) & (t<hmax))),over=over);
}

func hist(t,over=,color=, pasZero=, marks=)
/* DOCUMENT hist(tab,over=,color=)
        displays the histogram of array 'tab'.
        See also histo, or histogramme.

        Red  = median value
        Blue = average value
 */
{
  local t;

  if (pasZero) {
    nn = where(t!=0);
    if (is_array(nn)) t = t(nn);
  }
  
  limits,square=0;
  if( is_void(over) ) fma;
  logxy,0,0;
  t = double(t);
  if (min(t)!=max(t)) {
    y = histogram(int(1.5+254.*(t-min(t))/(max(t)-min(t))));
    x = span(min(t),max(t),dimsof(y)(2));
    tab = where(y>(max(y)/2));
    plg, y,x,color=color, marks=marks;
  }
//aaa = avg(t)
//mmm = median(t(*))
//plg,[0,max(y)],[mmm,mmm],marks=0,color="red"
//plg,[0,max(y)],[aaa,aaa],marks=0,color="blue"
// x(max(tab))-x(min(tab));
}

func gplh(y,x,color=,marks=,width=,type=)
{
  yy=xx=array(0.0, 2*numberof(y));
  yy(1:-1:2) = yy(2:0:2) = y;
  xx(2:-2:2) = xx(3:-1:2) = (x(2:-1))(pcen);
  xx(1) = x(1);
  xx(0) = x(0);
  plg, yy, xx, color=color, marks=marks,width=width,type=type;
}



 func colorBar(cmin, cmax)
/* DOCUMENT colorBar
            colorBar, cmin, cmax
     Draw a color bar to the right of the plot.  If CMIN and CMAX
     are specified, label the top and bottom of the bar with those
     numbers.

     EXAMPLE d'utilisation:
      pli,[span(.25,.75,200)] ;
      colorbar, .25,.75 ;
      
 */
{
  plsys, 0;
  pli, span(0,1,200)(-,), .625,.46,.67,.84, legend="";
  plg, [.46,.84,.84,.46],[.67,.67,.625,.625], closed=1,
    marks=0,color="fg",width=1,type=1,legend="";
  plsys, 1;  /* assumes there is only one coordinate system */
  if (!is_void(cmin)) {
    plt, pr1(cmin), .6475,.46, justify="CT";
    plt, pr1(cmax), .6475,.84, justify="CB";
  }
}

func clr(void)
{
  fma;
  limits;
}
func fullVisuMap( Caa )
/* DOCUMENT fullVisuMap( tmp ,ccc=, obs=)
   to compute maps from covariance matrices
   ccc = 1 means Cmaa computation
   ccc = 0 or void means Caa computation.
   
   SEE ALSO: visuMap
 */
  
{
  
  tmpx = dimsof(Caa)(2);
  tmpy = dimsof(Caa)(3);

  nssp = wfs(1).nLenslet;
  nbx = wfs(1).nValidSubAp;
  n = nssp*2-1;

  tmpx /= wfs(1).nValidMeas;
  tmpy /= wfs(1).nValidMeas;

  map = array(0.,(nssp*2-1)*2*tmpx,(nssp*2-1) *2*tmpy);
    
  for(i=1;i<=tmpx;i++) {        // remplissage de la metamap Caa
    for(j=1;j<=tmpy;j++) {
      indmapi = (i-1)*n*2 + 1:i*n*2;
      indmapj = (j-1)*n*2+1:(j)*n*2;
      indmati = (i-1)*nbx*2 + 1:i*nbx*2;
      indmatj = (j-1)*nbx*2+1:(j)*nbx*2
      map(indmapi,indmapj) = visuMap( Caa(indmati,indmatj), j);
    }
  }

  return map;
}

func visuMap( Caa, icam)
{

  nssp = wfs(icam).nLenslet;
  nbx = wfs(icam).nValidSubAp;
  n = wfs(icam).nLenslet*2-1;
  
  map = array(0.,n*2,n*2);
  map(1:n,1:n) = getCovMap(Caa(1:nbx,1:nbx),nssp); //yy
  map(1:n,n+1:) = getCovMap(Caa(1:nbx,nbx+1:),nssp ); // xy
  map(n+1:,n+1:) = getCovMap(Caa(nbx+1:,nbx+1:),nssp); //
  map(n+1:,1:n) = getCovMap(Caa(nbx+1:,1:nbx),nssp);

  return map;
}
func getCovMap(Caa, nssp)
/* DOCUMENT
     
   SEE ALSO:
 */

{

  // creation du masque des ssp
  x = span(-1,1,nssp)(,-:1:nssp);
  y = transpose(x);
  r = sqrt(x^2+y^2);
  param = 1.07;
  msk = r<param & r>tel.obs;
  nn = where(msk);
  // indices de la carte 2D carree de taille nsspXnssp, ou il y a 'vraiment' des sous-pupilles


  // creation du tableau des decalages
  xx = indgen(nssp)(,-:1:nssp)(nn); // coord en x des ssp
  dx = xx(,-) - xx(-,); // Ecart sur x entre chaque sous pupille
  yy = indgen(nssp)(-:1:nssp,)(nn); // coord en y des ssp
  dy = yy(,-) - yy(-,); // Ecart sur y entre chaque sous pupille

  // transformation des decalages en indice de tableau
  dx += nssp; 
  dy += nssp; 

  // transformation d'un couple de decalages (dx,dy) en un indice du tableau 'Map'
  map = div = array(0.0, nssp*2-1, nssp*2-1);
  ind = dx(*)+(nssp*2-1)*(dy(*)-1);
  for(i=1;i<=numberof(ind);i++) {
    map(ind(i)) += Caa(*)(i);
    div(ind(i)) += 1;
  }
  div( where(div==0) ) = 1;  // pour eviter division par 0
  map /= div;

  return map;
}
func plz(zer, color=, type=, marks=)
/* DOCUMENT  plz, varzer, color=, type=, marks=

     plz, array_varzer, color="black", type=2;
     will plot an array of variance on the zernike, in black with type 2 (dots)
     
     plz, 32.4
     will plot the Noll's spectrum on 35 Zernike, with the amplitude coefficient 32.4
   
   SEE ALSO:
 */
{
  noll = [0.448879,0.448879,0.0232179,0.0232179,0.0232179,0.00619143,0.00619143,
          0.00619143,0.00619143,0.00245392,0.00245392,0.00245392,0.00245392,0.00245392,
          0.00119041,0.00119041,0.00119041,0.00119041,0.00119041,0.00119041,0.000655102,
          0.000655102,0.000655102,0.000655102,0.000655102,0.000655102,0.000655102,
          0.000393378,0.000393378,0.000393378,0.000393378,0.000393378,0.000393378,
          0.000393378,0.000393378];
  
  if(isscalar(zer))
    zer = sqrt(noll)*zer;
  
  // number of the last zernike
  J = 1+dimsof(zer)(2);
  n = int( (-1.+sqrt(8*(J-1)+1))/2.); // radial order
  J2=1;
  for(i=1;i<=n;i++) {
    if(i>1) {
      plg, type=3, zer(J2-1:J2), [J2,J2+1], color=color, marks=0;
    }
    J1 = J2+1;
    J2 = (i+1)*(i+2)/2;
    J2 = min(J2,J);
    plg, zer(J1-1:J2-1), indgen(J1:J2), color=color, type=type, marks=0;
    if( marks!=0 )
      plmk, zer(J1-1:J2-1), indgen(J1:J2), color=color, marker=4, msize=0.3;
  }
}

func sortLabel(label_array)
/*DOCUMENT
  A tomoLabel is "MOAO4L3N", "SCAO","GLAO1L3N"... and the function returns a array of all the tomoLabel
  without redondancy

  ["MOAO4L3N","MOAO4L3N","SCAO","MOAO4L3N","GLAO4L3N","SCAO"] gives ["GLAO4L3N","MOAO4L3N","SCAO"].

  All the tomoLabel are in res_all.slopetl.label
 */
{

  //sort tomoLabel with the first letter
  type = label_array(sort(label_array));
  //defines first met label
  wtype = where(type==type(1));
  aa = "";
  //loop to retrieve the same tomoLabel
  while(1){
    aa = grow(aa,type(wtype(0)));
    if(wtype(0) != numberof(type)) wtype = where(type==type(wtype(0)+1));
    else break;
  }
  aa = aa(2:);

  return aa;
}

/*
 ____  _                       
/ ___|| | ___  _ __   ___  ___ 
\___ \| |/ _ \| '_ \ / _ \/ __|
 ___) | | (_) | |_) |  __/\__ \
|____/|_|\___/| .__/ \___||___/
              |_|              
               _                  _   _             
 _ __ ___  ___| |_ ___  _ __ __ _| |_(_) ___  _ __  
| '__/ _ \/ __| __/ _ \| '__/ _` | __| |/ _ \| '_ \ 
| | |  __/\__ \ || (_) | | | (_| | |_| | (_) | | | |
|_|  \___||___/\__\___/|_|  \__,_|\__|_|\___/|_| |_|
                                                    
*/

func determinesSlopesdisFromSlopestl(timetl)
/*DOCUMENT slopesdis = determinesSlopesdisFromSlopestl(timetl)

  Computes what the TS would see in open-loop from slopestl given by
  the local time time tl, voltages and MI. Returns the slopesdis from
  which we have remove the two first and last temporal samples. The
  data are reconstructed in pixel
 */
{



  slopestl = returnSlopestl(timetl,pathtl,arcsec=1);
  symTS = readFitsKeyArray(pathtl,"WFSSYM")(0);
  NslopesTS = readFitsKeyArray(pathtl,"WFSNSLO")(0);
  //loads voltages
  suffvolts = readFitsKey(pathtl,"VOLTFILE");
  if(suffvolts == "VOLTFILE not found")
    return 0;
  
  ptr_volts = restorefits("voltstl",suffvolts,pathvolts);
  if(dimsof(*ptr_volts(1))(1) == 0)
    return 0;
  //manages voltages
  voltsDM = *ptr_volts(1);//54xNframes matrix in voltages

  //loads Interaction Matrix
  suffMI = readFitsKey(pathtl,"MI");
  if(suffMI == "MI not found")
    return 0;
  
  MI = 3276.8*wfs(rtc.its).pixSize*restorefits("mi",suffMI,pathMI);//72x54 matrix in pixel/ADU
  
  //computes delay
  Fe = str2flt(readFitsKey(pathtl,"FREQ"));
  retard = 0.003*Fe + 1.05; //extrapolation
  fr= int(retard);
  coef = retard%1;
  
  slopesdis = slopestl;
  
  //convolution of the voltages
  vconv = coef*roll(voltsDM,[0,fr+1]) + (1.-coef)*roll(voltsDM,[0,fr]);
  //estimation of the corrected raw slopes
  rawslopes = MI(,+)*vconv(+,);
  //flip of the WFS
  slopesv = mirror_SH7(rawslopes, symTS);
  //adds to engaged slopes
  slopesdis(slrange(rtc.its),) -= slopesv;
  slopesdis = slopesdis(,3:-2);


  /*
  slopesdis = *rtc.slopes_res;
  volts = *rtc.volts;//in volts
  mi    = *rtc.mi; //in arcsec/volts

  frameDelay = rtc.frameDelay;
  fr= int(frameDelay);
  coef = frameDelay%1;

  //adds to engaged slopes
  if(rtc.obsMode != "LTAO"){
    //convolution of the voltages
    vconv = coef*roll(volts,[0,fr+1]) + (1.-coef)*roll(volts,[0,fr]);
    //estimation of the corrected raw slopes
    rawslopes = mi(,+)*vconv(+,); // in arcsec
    //flip of the WFS
    slopesv = mirror_SH7(rawslopes, wfs(rtc.its).sym);

    slopesdis(slrange(rtc.its),) -= slopesv;
    error;
  }else{
    //convolution of the voltages: to be determined for LTAO
    //vconv = coef*roll(volts,[0,*rtc.frameDelay+1]) + (1.-coef)*roll(volts,[0,*rtc.frameDelay]);
    //estimation of the corrected raw slopes
    rawslopes = mi(,+)*vconv(+,);
    //flip of the WFS
    slopesv = mirror_SH7(rawslopes, wfs(rtc.its).sym);
    
    slopesdis -= slopesv;
  }
  slopesdis = slopesdis(,3:-2);
  */
  return slopesdis;

}
func returnSlopestl(date,&path_data,fake=,arcsec =)
/*DOCUMENT slopestl = returnSlopestl(date,arcsec = )

  Returns slopestl from date in putting in the same format of the
  datatomo. If arcsec=1, returns the data in arcsec.
  
 */
{
  
  tmp = restorefits("slopestl",date,path_data,fake=fake);
  if(!fake){
    nwfs = numberof(tmp);
    slopestl = array(0., wfs(rtc.its).nValidMeas, nwfs, dimsof(*tmp(1))(4));
  
    for(i=1;i<=nwfs;i++){
      slopestl(,i,) = (*tmp(i))(*,);
    }
    for(i=1;i<=nwfs;i++) {
      slopestl(,i,) = mirror_SH7( slopestl(,i,) , wfs(i).sym );
      // remove WFS orientations for slopesdis
    }
    slopes_res = slopestl(*,);

    if(arcsec)
      for(i=1;i<=nwfs;i++) {
        slopes_res(slrange(i),) *= wfs(i).pixSize;
      }
    
    return slopes_res;
  }
    
}
func returnSlopesdis(date,&path_data,arcsec = )
/*DOCUMENT slopesdis = returnSlopesdis(date,arcsec = )

  Returns slopesdis from date in putting in the same format of the
  datatomo.  If arcsec=1, returns the data in arcsec.

 */
{
  
  tmp = restorefits("slopesdis",date,path_data);
  nwfs = numberof(tmp);
  slopesdis = array(0., wfs(rtc.its).nValidMeas, nwfs, dimsof(*tmp(1))(4));
  
  for(i=1;i<=nwfs;i++) slopesdis(,i,) = (*tmp(i))(*,);
  for(i=1;i<=nwfs;i++) {
    slopesdis(,i,) = mirror_SH7( slopesdis(,i,) , wfs(i).sym );// remove WFS orientations for slopesdis
  }
  slopes_dis = slopesdis(*,);

  if(arcsec)
    for(i=1;i<=nwfs;i++) {
      slopes_dis(slrange(i),) *= wfs(i).pixSize;
  }

  return slopes_dis;
    
}

func mirror_SH7(s_orig, sym, t=, inverse=)
/* DOCUMENT smir = mirror_SH7(s_orig, sym, t=, inverse=)

   SEE ALSO: mirror_SH_data
 */
{
  local inverse, t;
  return mirror_SH_data(s_orig,7,1.0,0.1,sym, t=t, inverse=inverse);
}
func mirror_SH_data(s_orig, n, radius, obs, sym, t=, inverse=)
/* DOCUMENT smir = mirror_SH_data(s_orig, nssp, ext_rad, obst_rad, sym, t=, inverse=)

   <s_orig> is an array of slopes, or an interaction matrix 
   <nssp> is the number of subapertures (7)
   <ext_rad> is the pupil external radius (1.07)
   <obst_rad> is the radius of central obstruction
   <sym> is the symmetry number, between 0 to 7, i.e.
         0 = nothing
         1 = symmetry around a vertical axis
         2 = symm around an horiz axis
         4 = symm around axis x=y (transpose)
   For several symmetries at once, just add.

   <t> allows to work on data that were transposed, i.e. where the dim
       of #subap is placed in 2nd position, like in interaction matrices.
   <inverse>=1 will apply the reciprocal transformation 
   
   SEE ALSO:
 */
{
  local inverse, t;
  sym%=8;
  if (!is_void(inverse) & inverse==1) {
    if (sym==5) sym=6;
    else if (sym==6) sym=5;
  }
  if( is_void(t) )
    t = 0;
  if(t==0)
    s = s_orig;
  else
    s = transpose(s_orig);
  valid = getValidSubapArray( n, radius, obs );
  nn = where(valid);
  msk = array(0,n,n);
  msk(nn) = indgen(1:numberof(nn));
  sm = s;
  xs=ys=1;
  if( (sym&1)==1 ) {
    msk = msk(::-1,);
    xs=-1;
  }
  if( (sym&2)==2 )  {
    msk = msk(,::-1);
    ys=-1;
  }
  if( (sym&4)==4 )  {
    msk = transpose(msk);
    ind = msk(nn);
    nssp = numberof(ind);
    sm(1:nssp,..) = s(ind+nssp,..) * ys;
    sm(nssp+1:,..) = s(ind,..) * xs;
  } else {
    ind = msk(nn);
    nssp = numberof(ind);
    sm(1:nssp,..) = s(ind,..) * xs;
    sm(nssp+1:,..) = s(ind+nssp,..) * ys;
  }
  
  if(t==0)
    return sm;
  else
    return transpose(sm);
}

/*
 ___           _           
|_ _|_ __   __| | _____  __
 | || '_ \ / _` |/ _ \ \/ /
 | || | | | (_| |  __/>  < 
|___|_| |_|\__,_|\___/_/\_\
                           
                                                             _   
 _ __ ___   __ _ _ __   __ _  __ _  ___ _ __ ___   ___ _ __ | |_ 
| '_ ` _ \ / _` | '_ \ / _` |/ _` |/ _ \ '_ ` _ \ / _ \ '_ \| __|
| | | | | | (_| | | | | (_| | (_| |  __/ | | | | |  __/ | | | |_ 
|_| |_| |_|\__,_|_| |_|\__,_|\__, |\___|_| |_| |_|\___|_| |_|\__|
                             |___/                               
*/

func slindex( icam )
/* DOCUMENT 

   Returns the index of the first slope of WFS number <icam> in an
   array where all the slopes are put one after the other also known
   as 'tomo format' with
   X-slopes1 Y-slopes1 X-slopes2 Y-slopes2 X-slopes3 Y-slopes3 ...
   When called as slindex( -icam ), the function returns the *last*
   slope of WFS number <icam>.

   slindex(1) = 1
   slindex(-1) = 72
   slindex(2) = 73
   slindex(-2) = 144
   slindex(3) = 145
   ...
   
   SEE ALSO:
 */
{
  if( icam<0 )
    return slindex(abs(icam)+1) - 1;
  n0 = 1;
  for(i=1;i<=icam-1;i++)
    n0 += wfs(i).nValidMeas;
  return n0;
}

func slrange(icam)
/* DOCUMENT 
   See tlsindex()
   
   SEE ALSO:
 */
{
  return slindex(icam):slindex(-icam);
}

func norange( x, n )
/* DOCUMENT norange( x )
            norange( x, n )
   Returns a range (whenever possible) or an array of indexes
   <x> is either
      - an array of WFS numbers (i.e. [1,3,rtc.its])
      - a scalar, i.e. rtc.its
      - a range, i.e. 1:rtc.nbLgs
      
   The function will return the list of index, between 1 and rtc.Nslopes, that are NOT
   within the range of the listed WFS numbers.
   
   SEE ALSO:
 */
{
  if(is_void(n)) n=sum(rtc.nSlopes);
  if(is_void(x))
    return 1:n;
  if( is_range(x) )
    x=indgen(x);
  if(isscalar(x)) x=[x];
  all = indgen(n);
  for(i=1; i<=numberof(x); i++) {
    k = x(i);
    all( slrange(k) )=0;
  }
  nn = where(all);
  if( is_array(nn) ) {
    if(allof(nn(dif)==1))
      return nn(1):nn(0);
    else
      return nn;
  } else
    return [];
}

func tsubrange(icam)
/* DOCUMENT 
   Returns the range of subapertures (not slopes!) for WFS number <icam> in an
   array where all the subapertures are put one after the other,
   X-subap1 Y-subap1 X-subap2 Y-subap2 X-subap3 Y-subap3 ...

   SEE ALSO:
 */
{
  n0 = 0;
  if(icam>1)
    n0 = sum(wfs(1:icam-1).nValidSubap);
  return 1+n0:n0+wfs(icam).nValidSubap;
}




func tslindex(icam)
/* DOCUMENT 
   Returns the index of the first slope of WFS number <icam> in an
   array where all the slopes are put one after the other also known
   as 'tomo format' with
   X-slopes1 Y-slopes1 X-slopes2 Y-slopes2 X-slopes3 Y-slopes3 ...
   When called as slindex( -icam ), the function returns the *last*
   slope of WFS number <icam>.

   tslindex(1) = 1
   tslindex(-1) = 72
   tslindex(2) = 73
   tslindex(-2) = 144
   tslindex(3) = 145
   ...
   
   SEE ALSO:
 */
{
  if( icam<0 )
    return tslindex(abs(icam)+1) - 1;
  n0 = 1;
  for(i=1;i<=icam-1;i++)
    n0 += wfs(i).nValidMeas;
  return n0;
}


func tslrange(icam)
/* DOCUMENT 
   See tlsindex()
   
   SEE ALSO:
 */
{
  return tslindex(icam):tslindex(-icam);
}
func others(i,n)
/* DOCUMENT o = others(i) or o = others(i,n)
   Returns all the elements of indgen(n) that are not in array <i>.
   By default, <n> = rtc.nbWfs.
   SEE ALSO:
 */
{
  if(is_void(n))
    n = rtc.nWfs;
  if( is_void(i) ) 
    return indgen(1:n);
  nn = where(i<=n);
  if( is_array(nn) )
    i = i(nn);
  else
    return indgen(1:n);
  nn = where(i>0);
  if( is_array(nn) )
    i = i(nn);
  else
    return indgen(1:n);
  u = indgen(n);
  u(i) = 0;
  nn = where(u);
  if( is_array(nn) ) return u(nn);
  else return [];
}
/*
 _____ _                
|_   _(_)_ __ ___   ___ 
  | | | | '_ ` _ \ / _ \
  | | | | | | | | |  __/
  |_| |_|_| |_| |_|\___|
                        
                                                             _   
 _ __ ___   __ _ _ __   __ _  __ _  ___ _ __ ___   ___ _ __ | |_ 
| '_ ` _ \ / _` | '_ \ / _` |/ _` |/ _ \ '_ ` _ \ / _ \ '_ \| __|
| | | | | | (_| | | | | (_| | (_| |  __/ | | | | |  __/ | | | |_ 
|_| |_| |_|\__,_|_| |_|\__,_|\__, |\___|_| |_| |_|\___|_| |_|\__|
                             |___/                               
*/

func time2str( x, sep )
/* DOCUMENT time = time2str( strtime, separ )
   Examples :
     > time2str( -1.234 )
     "22h45m58s"
     > time2str( -1.234 ,[":",":",""])
     "22:45:58"
     > time2str( -1.234 ,["heures ","minutes ","sec "])
     "22heures 45minutes 58sec "
     >

   SEE ALSO:
 */
{
  if( is_void(sep) ) sep=["h","m","s"];
  nbofx = numberof(x);
  if(nbofx==1)
    x = array(x,1);
  res = array(string, 1, nbofx);
  for(i=1;i<=nbofx;i++){
    y = x(i)%24;
    if(y<0 ) y+=24.;
    y *= 3600;         // translation into seconds
    y = long(y+0.5);   // round to nearest integer second
    h = y/3600;        // hours !
    y -= h*3600;       // remaining seconds
    m = y/60;
    s = y - m*60;
    h = h%24;
    res(,i) = swrite(format="%02d%s%02d%s%02d%s",int(h),sep(1),int(m),sep(2),int(s),sep(3));
  }
  if(nbofx==1)
    return res(1);
  else
    return res;
}

func str2time(str)
/* DOCUMENT h = str2time(str)
     Converts  "20h32m12s" into 20.53667
     Converts [["00h30m00s","00h45m00s"],["00h30m00s","00h30m00s"]] into [[0.5,0.75],[0.5,0.5]]
     
   SEE ALSO:
 */
{
  if( dimsof(str)(1)==0 ) {
    h=m=s=0;
    n = strfind(":", str, n=10);
    if(is_array(where(strpart(str,n)==":")))
      sread, str, format="%d:%d:%d ",h,m,s; // SLODAR-like hour format
    else
      sread, str, format="%dh%dm%ds ",h,m,s; // STYC-like hour format
    hh = h+m/60.+s/3600.;
    if(hh>12) hh-=24.;
    return hh;
  }
  hh = array(0.0, dimsof(str));   // array with same shape as str
  tmp = hh(*);
  tmpstr = str(*);
  for(i=1;i<=numberof(tmp);i++) {
    tmp(i) = str2time(tmpstr(i));
  }
  hh(*) = tmp;
  return hh;
}

func extractTime(pathtl)
{

  date = extractDate(pathtl);

  return strpart(date,-8:);
}

func superExtractDate(str)
{
  ns = numberof(str);
  str_ssj = array(string,ns);
  for(i=1;i<=ns;i++){
    str_ssj(i) = extractDate(str(i));
  }
  return str_ssj;
}
func extractDate( fullName )
/* DOCUMENT datetime = extractDate( fullName )
     Allows to extract the date and time fields from an absolute filename
     like "/Users/gendron/truc/Calibs_002/mi_001_2010-11-12T14h03m08_maMatrice.fits".

     Works on 'old' format 2010-11-12_14h03m08s and new FITS format 2010-11-12T14:03:08
     
   SEE ALSO:
 */
{
  a = decoupe(fullName)(0);   // filename alone without full path
  datetime = findPatternInString( a, "....-..-.._..h..m..s" );
  if(datetime=="")
    return findPatternInString( a, "....-..-..T..:..:.." );
  else
    return datetime;
}

func decoupe(str, delim)
/* DOCUMENT tab = decoupe(str, delim)
     decoupe("/users/canary/bin/")       produces ["users","canary","bin"]
     decoupe("/users/canary/bin/",'/');  produces ["users","canary","bin"]

     decoupe("/users/data/data_10-3-5T12:24:33_cui.fits",'_') produces
       ["/users/data/data","10-3-5T12:24:33","cui.fits"]
     
   SEE ALSO:
*/
{
  if( is_void(delim) )
    delim='/';
  a = strchar(str);
  nn = where(a==delim);
  if(is_array(nn))
    a(nn)=0;
  else
    return str;
  tabstr = strchar(a);
  nn = where(tabstr!=string(0));
  if(is_array(nn))
    return tabstr(nn);
  else
    return "";
}

func findPatternInString( str, pattern )
/* DOCUMENT
   > findPatternInString( "eric gendron", "g.n..o" )
   "gendro"
   > 
   > findPatternInString("datatomo_2011-07-03_19h54m52s_.fits",  "....-..-.._..h..m..s")
   "2011-07-03_19h54m52s"
   > 
     
   SEE ALSO:
 */
{
  u = strchar(pattern)(:-1);
  nu = numberof(u);
  s = strchar(str)(:-1);
  ns = numberof(s);
  for(i=1; i<=(ns-nu+1); i++) {
    v = s(i:i+nu-1);
    if(allof(v==u | u=='.')) return strchar(v);
  }
  return "";
}

func findTruc(truc,fic, margin=, after=, all=)
/* DOCUMENT findTruc(truc,fic, margin=, after=, all=)
            findTruc("offsetDM","23h50m03s", after=1)
            
   <margin> pour prendre une marge du fichier retrouve de "margin minutes".
   Si margin = 1 alors marge de 1 minute apres l'heure du fichier demande.

   <after> defaults to 0 (retrive file that *preceeds* requested time)
   <all> defaults to 0

   Si non specifie alors aucune marge n'est prise l'heure du fichier sera inferieure ou egale a celle specifiee
     
   SEE ALSO:
 */
{
  ficList = tmpDir + "list.txt";
  system, "ls -1tr "+dataDirRoot+"/"+truc+ "| cut -d '_' -f 3 > "+ficList;
  system_chmod,666,ficList;
  htime = rdfile(ficList);
  //error;
  if(!is_void(margin))
    margin=margin/60.;
  else
    margin=0.;
  if(is_void(htime)){
    swrite(format="ERROR CANNOT find %s. Exiting NOW!", truc);
    exit;
  }
  tmp = sort(str2time(htime)); // on remet ds l'ordre les fichiers
  time = str2time(htime)(tmp);
  if(!after)
    nn = where(time<str2time(fic)+margin);
  else
    nn = where(time>str2time(fic)+margin);
  
  if( is_array(nn) ) {
    if(!after){
      if(all)
        return time2str( time(nn) );
      else
        return time2str( time(nn)(0) );
    } else {
      if(all)
        return time2str( time(nn) );
      else
        return time2str( time(nn)(1) );
    }
  } else {
    return [];
  }
}

func checkRangeArray( tab, mini=, maxi=, type=, comment=, function=, dims=, nx=, ny=, nz= , file=, errorMessage=)
/* DOCUMENT  checkRangeArray( tab, mini=, maxi=, type=, comment=, function=, dims=, nx=, ny=, nz= )
   In the array tab, check that :
      - no one value is lower than mini
      - no one is greater than maxi
      - type is correct
      - number of dimensions dims is correct
      - first (nx) second (ny) and third (nz) dimensions are correct
    If not, display a message including "comment" and "function" values
      
   SEE ALSO:
 */
{
  tx = typeof( tab );
  errCode = 0;
  
  str = "Array";
  if( !is_void(comment) ) {
    str = comment;
  }

  buf = "";
  if( !is_void(function) ) {
    buf = "In function "+function;
    if (!is_void(file)) buf +=" ("+file+") : ";
    else buf+=" : ";
  }

  if (is_void(errorMessage)) errorMessage=0; 

  if( !is_void(type) ) {
    if( type!=tx ) {
      msg=swrite(format="%s%s is of type %s, while %s required.",buf,str,tx,type);
      errCode = 2;
    }
  }
  
  if( !is_void(mini) ) {
    nn = where(tab<mini);
    if( is_array(nn) ) {
      msg=write(format="%s%s contains values lower than %g.",buf,str,double(mini));
      errCode = 1;
    }
  }

  if( !is_void(maxi) ) {
    nn = where(tab>maxi);
    if( is_array(nn) ) {
      msg=swrite(format="%s%s contains values greater than %g.",buf,str,double(maxi));
      errCode = 1;
    }
  }
  
  if( !is_void(dims) ) {
    idim = dimsof(tab)(1);
    if( idim!=dims ) {
      msg=swrite(format="%s%s has %d dimensions (%d expected).",buf,str,idim,long(dims));
      errCode = 1;
    }
  }
  if (errCode) {
    if (errorMessage) return msg; else {write,msg+"\n"; return errCode;}
  }
  
  if( !is_void(nx) ) {
    idim = dimsof(tab);
    if (idim(1) < 1) {
      msg=swrite(format="%s%s has wrong dimension number : %d instead of at least 1 expected.",buf,str,idim(1));
      errCode = 1;
    }
    else {
      idim = idim(2);
      if( idim!=nx ) {
        msg=swrite(format="%s%s has wrong size : %d instead of %d expected.",buf,str,idim,long(nx));
        info, tab;
        errCode = 1;
      }
    }
  }
  if( !is_void(ny) ) {
    idim = dimsof(tab);
    if (idim(1) < 2) {
      msg=swrite(format="%s%s has wrong dimension number : %d instead of at least 2 expected.",buf,str,idim(1));
      errCode = 1;
    }
    else {
      idim = idim(3);
      if( idim!=ny ) {
        msg=swrite(format="%s%s has wrong size : %d instead of %d expected.",buf,str,idim,long(ny));
        info, tab;
        errCode = 1;
      }
    }
  }
  if( !is_void(nz) ) {
      idim = dimsof(tab);
    if (idim(1) < 3) {
      msg=swrite(format="%s%s has wrong dimension number : %d instead of at least 3 expected.",buf,str,idim(1));
      errCode = 1;
    }
    else {
      idim = idim(4);
      if( idim!=nz ) {
        msg=swrite(format="%s%s has wrong size : %d instead of %d expected.",buf,str,idim,long(nz));
        info, tab;
        errCode = 1;
      }
    }
  }
  if (errCode) {
    if (errorMessage) return msg; else {write,msg+"\n"; return errCode;}
  }
  return errCode;
}

func checkRange( x, mini=, maxi=, type=, comment=, function=, file=, errorMessage=)
/* DOCUMENT checkRange( x, mini=, maxi=, type=, comment=, function= )
   Check that x :
      - isn't lower than mini
      - isn't greater than maxi
      - type is correct
    If not, display a message including "comment", "function" and "file" values
      
   SEE ALSO:
 */
{
  errCode = 0;
  
  str = "Value";
  if( !is_void(comment) ) {
    str = comment;
  }

  buf = "";
  if( !is_void(function) ) {
    buf = "In function "+function;
    if (!is_void(file)) buf +=" ("+file+") : ";
    else buf+=" : ";
  }

  if (is_void(errorMessage)) errorMessage=0; 
  
  tx = typeof(x);
  if( tx=="complex" | tx=="struct_instance" | tx=="void" | tx=="range" | tx=="struct_definition" | tx=="function" | tx=="builtin") {
    msg = swrite(format="%s%s is of type %s : unsupported.",buf,str,tx);
    errCode = 4;
    if (errorMessage) return msg; else {write,msg+"\n"; return errCode;}
  }
  
  if( !is_void(type) ) {
    if( typeof(x)=="string" ) {
      if( type!="string" ) {
        msg = swrite(format="%s%s %s is of a wrong type (%s, when %s expected).",buf,str,x,typeof(x),type);
        errCode = 3;
        if (errorMessage) return msg; else {write,msg+"\n"; return errCode;}
      }
    }
    else {
      if( typeof(x)!=type ) {
        msg = swrite(format="%s%s %g is of a wrong type (%s, when %s expected).",buf,str,double(x),typeof(x),type);
        errCode = 2;
      }
    }
  }

  if( !is_void(mini) & errCode<3) {
    if( x<mini ) {
      msg = swrite(format="%s%s %g is lower than %g.",buf,str,double(x),double(mini));
      errCode = 1;
    }
  }

  if( !is_void(maxi) & errCode<3 ) {
    if( x>maxi ) {
      msg=swrite(format="%s%s %g is greater than %g.",buf,str,double(x),double(maxi));
      errCode = 1;
    }
  }
  if (errCode) {
    if (errorMessage) return msg; else write,msg+"\n";
  }
  return errCode;
}


func E2880(i)
/* DOCUMENT j = E2880(i)
     returns the first multiple of 2880 greater or equal than i
   SEE ALSO:
 */
{
  return (1l+(long(i)-1l)/2880l)*2880l;
}


/*
 _____                          _   
|  ___|__  _ __ _ __ ___   __ _| |_ 
| |_ / _ \| '__| '_ ` _ \ / _` | __|
|  _| (_) | |  | | | | | | (_| | |_ 
|_|  \___/|_|  |_| |_| |_|\__,_|\__|
                                    
                                                             _   
 _ __ ___   __ _ _ __   __ _  __ _  ___ _ __ ___   ___ _ __ | |_ 
| '_ ` _ \ / _` | '_ \ / _` |/ _` |/ _ \ '_ ` _ \ / _ \ '_ \| __|
| | | | | | (_| | | | | (_| | (_| |  __/ | | | | |  __/ | | | |_ 
|_| |_| |_|\__,_|_| |_|\__,_|\__, |\___|_| |_| |_|\___|_| |_|\__|
                             |___/                               
*/

func int2str(i)
/* DOCUMENT my_string = int2str(i)
     Converts an integer to a string.
   SEE ALSO:
 */
{
  if( checkRange(i,comment="integer value",function="int2str") )
    exit;
  return swrite(format="%d",long(i));
}

func var2str(xvalue)
/* DOCUMENT str = var2str(x)
     Converts an array of numbers (floats, double or long) into string array.
     Keeps a
   SEE ALSO:
 */
{
  tx=typeof(xvalue);
  if(tx=="string") return xvalue;
  value = [];
  if(tx=="char"|tx=="int"|tx=="long"|tx=="short") {
    value = swrite(format="%d",long(xvalue));
  }
  if(tx=="float"|tx=="double") {
    value = swrite(format="%g",double(xvalue));
  }
  return value;
}

func str2flt(s, noerr=)
/* DOCUMENT str2flt(s, noerr=)
     Converts a string like "3.14159" to a floating point value 3.14159.
     The function makes an error when attemping to convert strings
     like "toto".
     When optional parameter <noerr>=1 it returns 0.00 with no error.
     
   SEE ALSO:
 */
{
  x = array(0.0, dimsof(s));
  tmps = sum(s+" ");
  tmp = x(*);
  n = sread(tmps,format="%f ",tmp);
  if( n<1 ) {
    if( noerr==1 )
      return 0;
    else
      error, "Can't convert string \""+s+"\" into float. Exiting.";
  }
  x(*) = tmp;
  return x;
}

func str2int(s)
{
  x = array(0, dimsof(s));
  tmps = sum(s+" ");
  tmp = x(*)
  sread,tmps,format="%d ",tmp;
  x(*) = tmp;
  return x;
}

func anything2str(xvalue)
/* DOCUMENT str = anything2str(x)
     Converts an array of numbers (floats, double or long) into string array.
     
   SEE ALSO:
 */
{
  tx=typeof(xvalue);
  if(tx=="string") return xvalue;
  value = [];
  if(tx=="char"|tx=="int"|tx=="long"|tx=="short") {
    value = swrite(format="%d",long(xvalue));
  }
  if(tx=="float"|tx=="double") {
    value = swrite(format="%g",double(xvalue));
  }
  return value;
}

func notExist(a)
/* DOCUMENT notExist(a)
  return 1 if the file does not exist, else 0
*/
{
  fic = open( a, "rb", 1);
  if( fic ) {
    close, fic;
    return 0;
  }
  return 1;
}

func isscalar(x)
/* DOCUMENT isscalar(x)
     return 1 if x is a scalar, 0 else
   SEE ALSO:
 */
{
  if( is_void(x) )
    return 0;
  if( x==where() )
    return 0;
  if( is_array(x) & dimsof(x)(1)==0 & kindoftype(x)=="numeric" )
    return 1;
  else
    return 0;
}
func kindoftype(x)
/* DOCUMENT kindoftype(x)
     return the type of x : "numeric" or "not_numeric"
   SEE ALSO:
 */
{
  tx=typeof(x);
  if(tx=="char"|tx=="int"|tx=="long"|tx=="short"|tx=="float"|tx=="double") {
    return "numeric";
  }
  else {
    return "not_numeric";
  }
}

func file_exists(filename)
{
  fic = open(filename,"r",1); 
  if (fic) {
    close,fic;
    return 1;
  }
  else{
    return 0;
  }
}
func direxist(a)
{
  // definition of the file that will contain the list of the directory
  ficlist = tmpDir+"/lsList.txt";
  // rm things
  com_rm = "rm -f "+tmpDir+"dialogFile.txt "+ficlist+" > /dev/null 2>&1 ; ";
  system,com_rm + "ls -ld "+a+" > "+ficlist+" 2> /dev/null ; echo $? > "+tmpDir+"dialogFile.txt ; chmod 666 "+tmpDir+"dialogFile.txt > /dev/null 2>&1";
  // read the file to know if 'ls' got no error
  str = rdfile( tmpDir+"dialogFile.txt" )(0);
  if( str=="0" ) {
    // file exist. Now checking if it's really a directory
    fic = open(ficlist,"rb",1);
    if(fic==[]) return 0;
    car='-';
    _read,fic,0,car;
    close,fic;
    if(car=='d') return 1;
    return 0;
  } else {
    return 0;
  }
}
/*
               _                  _   _             
 _ __ ___  ___| |_ ___  _ __ __ _| |_(_) ___  _ __  
| '__/ _ \/ __| __/ _ \| '__/ _` | __| |/ _ \| '_ \ 
| | |  __/\__ \ || (_) | | | (_| | |_| | (_) | | | |
|_|  \___||___/\__\___/|_|  \__,_|\__|_|\___/|_| |_|
                                                    
*/

func listFile(path)
{

  // definition of the file that will contain the list of the directory
  ficlist = tmpDir+"lsList.txt";
  // first, rm of the file (in order that it's not confused with a previous call)
  system, "rm -f " + tmpDir + "dialogFile.txt " + ficlist + " > /dev/null 2>&1";
  // file listing, in reverse time order
  system, "ls -1tr " + path + "> " + ficlist + " 2> /dev/null ; echo $? > " + tmpDir + "dialogFile.txt";
  // put permissions for all the newly created files
  system_chmod,666,tmpDir+"dialogFile.txt";
  system_chmod,666,ficlist;
  // read the file to know if 'ls' got no error
  str = rdfile( tmpDir+"dialogFile.txt" )(0);
  if( str=="0" ) {
    listdir = rdfile(ficlist);
    return listdir;
  } else {
    return [];
  }
}

func listVersion( root, type, filePrefix, fileSuffix )
{
  if( is_void(fileSuffix) ) fileSuffix="";
  if( fileSuffix==" " ) fileSuffix="";    // to avoid a biiiig bug .....!!!
  // definition of the file that will contain the list of the directory
  ficlist = tmpDir+"/lsList.txt";
  // first, rm of the file (in order that it's not confused with a previous call)
  system, "rm -f "+tmpDir+"dialogFile.txt "+ficlist+" > /dev/null 2>&1";
  // file listing, in reverse time order
  system, "ls -1tr "+root+"/"+filePrefix+"/"+filePrefix+"*"+fileSuffix+"*."+type+" > "+ficlist+" 2> /dev/null ; echo $? > "+tmpDir+"dialogFile.txt";
  // put permissions for all the newly created files
  system_chmod,666,tmpDir+"dialogFile.txt";
  system_chmod,666,ficlist;
  // read the file to know if 'ls' got no error
  str = rdfile( tmpDir+"dialogFile.txt" )(0);
  if( str=="0" ) {
    // ls got not error ...
    listdir = rdfile(ficlist);
    return listdir;
  } else {
    // error in 'ls' ... no file ?
    return [];
  }
}

func system_mkdir( dirname )
{
  system, "mkdir "+dirname+" ; chmod 777 "+dirname+" > /dev/null 2>&1";
}

func system_chmod( x , filename )
{
  com = "chmod "+int2str(x)+" "+filename+" > /dev/null 2>&1";
  system,com;
}
func readFile(filename,type,fake=)
{
  local fake;
  extern tomo;
  if( strmatch(type,"fits",1) ) { // FITS type file !
    if( !fake ) {
      dat = READFITS(filename);
      if( dat==[] ) exit;
      return dat;
    }
    return 0;
  }


  if( strmatch(type,"ybin",1) ) { // ybin type file of struct variable called <tomo>
    if( !fake ) {
      dat = restoreTomoStructure( filename );
      if( dat==[] ) exit;
      return dat;
    }
    return 0;
  }

  if( strmatch(type,"txt",1) ) { // ybin type file of struct variable called <tomo>
    if( !fake ) {
      dat = scidar_readprofile(filename);
      // dat = readasciiMAC(filename);
      if( dat==[] ) exit;
      return dat;
    }
    return 0;
  }

  
  for(i=1;i<=100;i++)
    write," type is unknown !!!! can't read this !! ";
  error;
}

func restorefits( filePrefix, fileSuffix, &fullName, version=, fake=,verb= )
/* DOCUMENT restorefits( filePrefix, fileSuffix, fullName, version=, fake= )
     Restores a FITS file from the archive.
     Argument <fullName> is an output, it contains the full path name.
     
     Examples:
     a = restorefits("caa");
     a = restorefits("caa", "good");
     a = restorefits("caa", version=-2);
     a = restorefits("caa",,fullName);
     a = restorefits("caa",,fullName, version=i, fake=1);
     
   SEE ALSO:
 */
{
  local version, fake;
  extern dataDirRoot;
  // check that filePrefix is not void
  if( is_void(filePrefix) ) {
    helpkeyword;
    return [];
  }
  // if everything's ok, restore the data ....
  data = restoreContext( dataDirRoot, "fits", filePrefix, fileSuffix, fullName, version=version, simu=1, fake=fake );

  // check the data format wrt STYC version when slopestl or caa are required
  if( data!=[] )
    data = manageDataFormatWrtStycVersion(data, filePrefix, fullName, fileSuffix, version, fake);
  
  return data;
}

func restoreContext( root, type, filePrefix, fileSuffix, &fullName, version=, fake=, simu= )
/* DOCUMENT 
     restoreContext( root, type, filePrefix, fileSuffix, &fullName, version=, fake=, simu= )
     
   SEE ALSO:
 */
{
  local version, fake, simu;

  if( is_void(version) ) version=0;
  if( is_void(simu) ) simu=0;
  if( is_void(fake) ) fake=0;
  if( is_void(fileSuffix) ) fileSuffix="";
  
  // listing proper directory
  str = listVersion( root, type, filePrefix, fileSuffix );
  // <nfound> files have been found
  nfound = numberof( str );
  // check <nfound> is non-zero
  if( nfound ) {
    // check that number <version> can be retrieved
    if( (version>nfound) || (version<(1-nfound)) ) {
      fullName = swrite(format="ERROR : impossible to retrieve data version %d. Only %d versions available", version, nfound);
      write,format="%s\n",fullName;
      write,format="File searched was : %s\nwith suffix       : %s\n",filePrefix, fileSuffix;
      fullName = "/error/"+fullName;
      return [];
    }
    fullName=str(version);
    
    return readFile( fullName, type, fake=fake );
  } else {
    //write,format="ERROR : impossible to retrieve file %s.\n",filePrefix;
    //error;
      fullName = swrite(format="ERROR : impossible to retrieve data version %d. Only %d versions available", version, nfound);
      write,format="%s\n",fullName;
      write,format="File searched was : %s\nwith suffix       : %s\n",filePrefix, fileSuffix;
      fullName = "/error/"+fullName;
      return [];
  
  }
}
func savefits( data, filePrefix, fileSuffix )
/* DOCUMENT fullname = savefits( data, filePrefix, fileSuffix )
     
   SEE ALSO:
 */
{
  extern dataDirRoot;
  file = saveContext( dataDirRoot, data, "fits", filePrefix, fileSuffix, fake=0 );
  return file;
}
func saveContext( root, data, type, filePrefix, fileSuffix, fake= )
/* DOCUMENT saveContext( root, data, type, filePrefix, fileSuffix, fake= )
     
   SEE ALSO:
 */
{
  local fake, over;
  if( fileSuffix==[] ) fileSuffix="";
  if( is_void(data) ) {
    write,format="ERROR : saveContext() has detected the data to be saved (type=%s, prefix=%s, suffix=%s) in just VOID !!!\n",type,filePrefix,fileSuffix;
    write,"ERROR : saveContext() will not save anything and just return and empty string as pathname.";
    return "";
  }
  if( typeof(filePrefix)!="string" ) {
    helpkeyword;
    error;
  }
  dirsave = root+"/"+filePrefix;
  if( !direxist(dirsave) ) {
    // rm things
    com_rm = "rm -f "+tmpDir+"dialogFile.txt > /dev/null 2>&1 ; ";
    system, com_rm+"mkdir "+dirsave+" > /dev/null 2>&1 ; echo $? > "+tmpDir+"dialogFile.txt";
    // read the file to know if 'ls' got no error
    str = rdfile( tmpDir+"dialogFile.txt" )(0);
    if( str!="0" ) {
      error;
    }
  }

  fullName = writeFile(dirsave,data,type,filePrefix,fileSuffix,fake=fake);
  if(fake!=1)
    lastData, root, filePrefix, fullName, type;
  
  return fullName;
}
func thisIsTheLastFile( path )
{
  hiddenfile = dataDir+"/"+".thisIsTheLastFile";
  system, "echo "+path+" > "+hiddenfile;
  system_chmod,666,hiddenfile;
}
func addTheLastFile( filename )
{
  listOfAllFilesWritten = dataDirRoot+"/"+"listOfAllFilesWritten";
  system, "echo "+filename+" >> "+listOfAllFilesWritten;
  system_chmod,666,listOfAllFilesWritten;
}

func littlewrite( format, filename )
{
  // This is the function that writes all the filenames on the yorick screen ...
  // It's been indentified as a separate function, in order to be able to
  // disable it
  if( 1 )
    write,format=format, filename;
}

func writeFile(path,data,type,filePrefix,fileSuffix,fake=)
{
  local fake;
  extern tomo;
  extern KEY;
  
  // if no fileSuffix is given
  if( is_void(fileSuffix) )
    fileSuffix = "";
  // remove all spaces and other strange characters from string such as  [ ] , : /
  fileSuffix = nospc(fileSuffix);
  // filename = path+"/"+filePrefix+"_"+version+"_"+makeMyDay()+"_"+fileSuffix;  // ca, c'etait avant !!!!...
  filename = path+"/"+filePrefix+"_"+makeMyDay()+"_"+fileSuffix;

  thisIsTheLastFile, path;
    
  if( strmatch(type,"fits",1) ) {          // FITS type file !
    filename += ".fits";
    if( !fake && NO_FILE_WRITING!=1 ) {
      key = chooseKeywords(filePrefix);
      WRITEFITS, filename, data, key(1,), key(2,);
      system,"rm -f "+filename+"L > /dev/null 2>&1 ";
      littlewrite, "Writing :   \n%s\n", filename;
      run = decoupe(dataDirRoot,'/')(0);
      sss = decoupe(filename,'/');
      ind=where(sss==run)(1);
      filenameRel = strpart(makePath(sss(ind+1:)), 2:);
      addTheLastFile, filenameRel;
    }
    return filename;
  }
  
  if( strmatch(type,"txt",1) ) { // TXT file !
    filename += ".txt";
    if( !fake ) {
      asciiDump,filename,data;
      littlewrite, "Writing :   \n%s\n", filename;
    }
    return filename;
  }


  if( strmatch(type,"ybin",1) ) { // will save the tomo structure (or any other structure, by the way..)
    filename += ".ybin";
    if( !fake ) {
      //archiveTomoStructure, filename, data; @@@@@@@@@@@@@@@@@@@@@@@@@
      littlewrite, "Writing :   \n%s\n", filename;
    }
    return filename;
  }

  
  for(i=1;i<=100;i++)
    write,"AAAAAAAAAAAAAAAAAAAAAAAAAAAA";
  error;
}

func lastData(root, filePrefix, fullName, type)
/* DOCUMENT lastData(root, filePrefix, fullName, type)
     
   SEE ALSO:
 */
{
  if( !direxist(root+"/last") ) {   // Directory does NOT exist !!
    system_mkdir, root+"/last";
  }

  //  system,"cp -f "+fullName+" "+root+"/last/"+filePrefix+"."+type;
  system,"ln -fs "+fullName+" "+root+"/last/"+filePrefix+"."+type+" > /dev/null 2>&1 ";
}
func makePath( listDir )
/* DOCUMENT p = makePath( listDir )
   
     > makePath(["Users","gendron","tmp","truc","machin"])
     "/Users/gendron/tmp/truc/machin"
     >
     
   SEE ALSO:
 */
{
  p = "";
  n = numberof(listDir);
  for(i=1;i<=n;i++) p += "/"+listDir(i);
  return p;
}




/*

 ____        _          _____                          _   
|  _ \  __ _| |_ __ _  |  ___|__  _ __ _ __ ___   __ _| |_ 
| | | |/ _` | __/ _` | | |_ / _ \| '__| '_ ` _ \ / _` | __|
| |_| | (_| | || (_| | |  _| (_) | |  | | | | | | (_| | |_ 
|____/ \__,_|\__\__,_| |_|  \___/|_|  |_| |_| |_|\__,_|\__|
                                                           

               _   
__      ___ __| |_ 
\ \ /\ / / '__| __|
 \ V  V /| |  | |_ 
  \_/\_/ |_|   \__|
                   

 ____  _               __     __            _             
/ ___|| |_ _   _  ___  \ \   / /__ _ __ ___(_) ___  _ __  
\___ \| __| | | |/ __|  \ \ / / _ \ '__/ __| |/ _ \| '_ \ 
 ___) | |_| |_| | (__    \ V /  __/ |  \__ \ | (_) | | | |
|____/ \__|\__, |\___|    \_/ \___|_|  |___/_|\___/|_| |_|
           |___/                                          



 */

func manageDataFormatWrtStycVersion(data,filePrefix,fullName,fileSuffix,version,fake)
/* DOCUMENT data = manageDataFormatWrtStycVersion(data,filePrefix,fullName,fileSuffix,version,fake)

   Changes the format of the data read by styc from a FITS file, in
   order to be compatible between the different Styc versions.

   Styc version 2.0: Data need to be rearranged if they are slopestl
   or caa/cmaa written by Styc 1.0.
  
   SEE ALSO:
 */
{
  if( fake==1 ) return data;
  // data need to be rearranged if they are slopestl or caa/cmaa
  // because format has been changed between styc versions 1.0 and 2.0
  if( filePrefix=="slopestl" || filePrefix=="caa" || filePrefix=="caafit" || filePrefix=="abstats" || filePrefix=="tomoparam" ) {
    // reading styc version in the file header
    stycVersion = getStycVersion(fullName);
    if( stycVersion<2.0 ) {
      if( filePrefix=="slopestl" )
        return rearrangeSlopestl(data);
      if( filePrefix=="caa" ) {
        // In Styc 1.0 the covariance matrix was saved in 2 parts, in 2 different files : Caa, and Cmaa.
        cmaa = restoreContext( dataDirRoot, "fits", "cmaa", fileSuffix, fullNameCaa, version=version, simu=1, fake=fake );
        return rearrangeCaa(data, cmaa);
      }
      if( filePrefix=="caafit" ) {
        // In Styc 1.0 the covariance matrix was saved in 2 parts, in 2 different files : Caa, and Cmaa.
        cmaa = restoreContext( dataDirRoot, "fits", "cmaafit", fileSuffix, fullNameCaa, version=version, simu=1, fake=fake );
        return rearrangeCaa(data, cmaa);
      }
      if( filePrefix=="abstats" ) {
        return rearrangeAbstats(data);
      }
      if( filePrefix=="tomoparam" ) {
        return rearrangeTomoparam(data);
      }
    }
  }
  return data;
}
func rearrangeTomoparam(data)
{
  nbcam = 4;
  // ............. what follows is a copy of the function unpackcoeffs_fab() ..............
  nbCamOffaxis = nbcam-1;
  pt = 1; // indice de demarrage
  nonoise=0; // pour detecter si il y a pas de noise dans les coeffs.
  l = 3;               // longueur du nbre d'elements a lire dans le tableau
  cste = data(pt:pt+l-1);
  pt+=l;

  l = 1;
  L0 = data(pt:pt+l-1);
  pt+=l;
  
  nbl = ((dimsof(data)(2) - (nbcam-1) - 1-(nbcam-1)*2 -2 - (nbcam-1))/2);
  if(nbl==0){
    nbl = (dimsof(data)(2) - (nbcam-1) - 1-(nbcam-1)*2 -2)/2;
    nonoise=1;
  }
  l = 2*nbl;
  Cn2h = data(pt:pt+l-1);        // tableau contenant [h1,r01,h2,r02......hn,r0n]
  pt+=l;
  r0 = Cn2h(2::2);
  altitu = Cn2h(1::2);

  l = (nbCamOffaxis+1)*2;
  posStar = data(pt:pt+l-1);       // tableau contenant les positions relatives des ASOs par rapport au premier (pris c reference).
  pt+=l;
  xposStar = posStar(1::2);
  yposStar = posStar(2::2);

  if(!nonoise){
    l = nbCamOffaxis; // Bon ok valeurs de bruit sur la diag pour l'instant.
    noise = data(pt:pt+l-1);
    noise = grow(noise,0.0)
  } else {
    noise = array(0., nbcam);
  }
  // .................. end of unpackcoeffs_fab() ..............

  
  // ........ transformation to STYC2 format .........
  //
  //   WARNING : units are wrong, they need to be rescaled
  //   <altitu> is normalized, from 0 to 1 for the highest layer
  //   <r0> is proportionnal to a 'real' r0, not to a phase variance
  //   <xposStar> and <yposStar> are expressed in subapertures units
  //   L0 and noise ... ????
  //
  write,"WARNING : data have been transformed into STYC2 format ..";
  write,"          Please check units !!!!";
  write,"<altitu> is normalized, from 0 to 1 for the highest layer";
  write,"<r0> is proportionnal to a 'real' r0, not to a phase variance";
  write,"<xposStar> and <yposStar> are expressed in subapertures units";
  write,"L0 and noise ... ????";
  ptr = [ &L0, &r0, &altitu, &noise, &cste, &xposStar, &yposStar, &array(1.0,nbcam) ];
  return ptr;
}
func rearrangeAbstats(data)
{
  n = dimsof(data)(3);   // nber of WFSs
  ptr = array(pointer,n);
  for(i=1; i<=n; i++) ptr(i) = &(foldxy(data(,i)));
  return ptr;
}
func rearrangeCaa(caa, cmaa)
{
  nx = numberof(caa(,1)) + numberof(cmaa(,1)); // total size when caa and cmaa will be concatenated
  fullCaa = array(0.0, nx, nx);
  fullCaa(norange(rtc.its), norange(rtc.its)) = caa;
  fullCaa(norange(rtc.its), slrange(rtc.its)) = transpose(cmaa);
  fullCaa(slrange(rtc.its), norange(rtc.its)) = cmaa;
  return fullCaa;
}
func rearrangeSlopestl(data)
/* DOCUMENT data = rearrangeSlopestl(data)
     Will rearrange slopestl from version styc1.0 to styc2.0.

     - Styc 1.0 was a format (36,2,4,1000) while Styc 2.0 is an array
       of 4 pointers on (36,2,1000) arrays.
     - Sometimes in styc1.0, the optigain procedure of matthieu saved
       some slopestl in a format (72,1000).
     
   SEE ALSO:
 */
{
  s = dimsof(data);
  if(s(1)==2) {
    // only one WFS in a format 72 ... (usually created by 'optigain', in 2010)
    return foldxy(data);
  } else if(s(1)==4) {
    // 4 WFSs, format (36,2,4,1000)
    n = s(4);  // number of wfs
    ptr = array(pointer,n);
    for(i=1; i<=n; i++) ptr(i) = &(data(,,i,));
    return ptr;
  } else {
    // 1 WFS, format (36,2,1000) or anything else
    return data;
  }
}
func getStycVersion(fullName)
/* DOCUMENT stycVersion = getStycVersion(fullName);
   
   Returns the floating-point number corresponding to the styc version
   that wrote the fits file.
   If the file has been written before 2011, then no STYC keyword is
   present in the header, and the number 1.0 is returned.
   
   SEE ALSO:
 */
{
  stycVersion = readFitsKey(fullName, "STYC");
  if( strmatch(stycVersion,"not found") )
    stycVersion="1.0";
  stycVersion = str2flt(stycVersion);
  return stycVersion;
}


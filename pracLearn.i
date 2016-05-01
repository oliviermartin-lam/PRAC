/*
 _____ _ _   _   _             
|  ___(_) |_| |_(_)_ __   __ _ 
| |_  | | __| __| | '_ \ / _` |
|  _| | | |_| |_| | | | | (_| |
|_|   |_|\__|\__|_|_| |_|\__, |
                         |___/ 
*/
func fitCovarianceMatrix(slopesArray,nl,Cnn=,ttr=,fitl0=,FitMethod=,fullHD=,tomores=,verb=,deriv=,misreg=)
/* DOCUMENT ptrcn2h = fitCovarianceMatrix(*data.slopes_dis,3,fitl0=1,FitMethod=1,fullHD=0,tomores=5000,verb=1,deriv=0)

   Performs a quick fitting from input arguments.
 */
{
  learn.nl = nl;
  if(ttr) learn.ttr = 1;
  //.....INITIAL GUESS.....//
  tmp = initGuessLearn(nl,fitl0=fitl0,fullHD=fullHD,tomores=tomores,misreg=misreg);
  pinit = tmp(,1);
  pind = tmp(,2);

  //.....Compute covariance matrix of turbulence with square ssp.....//
  slopesArray -= slopesArray(,avg);
  Cdd = slopesArray(,+)*slopesArray(,+)/dimsof(slopesArray)(0);

  //fit
  if(is_void(FitMethod)){
    FitMethod = 1;
  }
  if(is_void(Cnn))
    Cnn = *covMatrix.noise;

  ptrcn2h = doLearn(Cdd,FitMethod,pinit,pind,Cnn,verb=verb);
 
  return ptrcn2h;
}

func findAltitudeMax(posX,posY,&Bmax,useTS=,altituderes=,altLGS=)
/* DOCUMENT altmax = findAltitudeMax(posX,posY,useTS=,altsampling=)

  Returns the maximal altitude seen by Canary in Phase B config. Pos X
  and Pos Y are the tas vector positions in arcsec. If useTS!=0, the
  maximal altitude is computing in taking into account the base
  between the TS and the off-axis WFS. If altituderes ==1, returns the
  resolution in altitude we get from the shortest base we have
  
   Olivier Martin.
 */
{

  //write,format="%s\n","THIS FUNCTION DOESN T WORK IN PHASE C !!!";
  //Define the number of WFS
  n = numberof(posX)-1;
  if(useTS) n+=1;
  
  posX = posX(1:n);
  posY = posY(1:n);

  if(is_void(altLGS))
    altLGS = 21000.;
  if(is_struct(tomo))
    alt_lgs = tomo.wfs.alt_lgs;
  else
    alt_lgs = array(altLGS,n);
  //Bases array
  B = [];
  D = tel.diam;
  //Base between each stars
  for(i=1;i<=n-1;i++){
    if(useTS || i!=n+1){
      Btmp = sqrt((posX(i) - posX(i+1:))^2 + (posY(i) - posY(i+1:))^2);
      if(alt_lgs(i)){
        Blim = D*radian2arcsec/alt_lgs(i);
        Btmp*=(Btmp > Blim);//if altitude max > LGS altitude: but baseline to 0
      }
      B = grow(B,Btmp);
    }
  }
  B = B(where(B!=0))/206265.;
  
  if(anyof(B)){
    Bmax = max(B);Bmin=min(B);
    if(altituderes)
      return  D/Bmax/wfs(rtc.its).nLenslet;
    else 
      return D/Bmin;
  }

}

func initGuessLearn(nbl,fitl0=,fullHD=,tomores=,misreg=)
  /* DOCUMENT [ptr_initarray,pind] = initGuessLearn(nb_layers,posX,posY,fullHD=,tomores=,misreg=)

     Renvoie deux pointeurs pour les conditions initiales du Learn
     (ptr_initarray) et les parametres a fitter (pind).
     nb_layers est le nombre de couches turbulentes et posX,posY les
     coordonnees du TAS et du TS en arcsec. Les altitudes sont
     reparties eqalement entre 0 et l'altitude max determinee par la
     geometrie de l'asterisme (sans TS). S fullHD==1, l'altitude des
     couches ne sera pas fittee. Si tomores est precisee, les
     altitudes seront reparties entre 0 et tomores*(nb_layers-1). Le
     cn2h est egal a 0.1 pour toutes les couches sauf au so ou il vaut
     1.  Les sensibilites sont fixees a 1, L0 a 1000, tracking et
     bruit a 0.


     Olivier Martin.

     SEE ALSO: findAltitudeMax
 */
  
{
  altmax = findAltitudeMax(wfs.x,wfs.y,useTS=1);//find the maximal altitude seen by Canary without TS
  pinit = pind = array(pointer,8);

  //1. cnh profile
  cnh = 20.;
  cnhfit = 1.
  if(nbl > 1){
    cnh = grow(20.,array(1.,nbl-1));
    cnhfit = array(1.,nbl);
  }
  pinit(1) = &cnh;
  pind(1) = &cnhfit;

  //2. altitude
  if(nbl == 1){
    alt = 0;
  }else{
    alt = span(0,altmax,nbl);
  }
  if(tomores && nbl!=1){
    alt = span(0,(nbl-1)*tomores,nbl);
  }

  pinit(2) = &alt;
  fitalt=1;
  if(fullHD){
    fitalt = 0;
  }
  pind(2) = &(array(fitalt,nbl));

  //3. l0 profile
  l0h = atm.L0;
  if(nbl > 1 && fitl0 == 2){
    l0h = array(10.,nbl);
  }else if(fitl0 !=2){
    l0h = grow(10.,array(100.,nbl-1));
  }
  pinit(3) = &l0h;
  if(is_void(fitl0)| fitl0 == 0){
    pind(3) = &(array(0.,nbl));
  }else if(fitl0 == 1){
    pind(3) = &(grow(1,array(0.,nbl-1)));
  }else if(fitl0 == 2){
    pind(3) = &(array(1.,nbl));
  }
  
  //4. Tracking
  pinit(4) = &array(0.,3);
  if(learn.ttr == 0){
    pind(4) = &array(1,3);
  }else{pind(4) = &array(0,3);}

  amis = array(0.,rtc.nWfs);
  if(misreg == 1)
    amis(1:-1) = 1;
  
  //shift pupil
  pinit (5) = pinit (6) = &array(0.,rtc.nWfs);
  pind (5) = pind (6) = &amis;

  //magnification
  pinit (7) =  &array(1.,rtc.nWfs);
  pind  (7) =  &amis;

  //rotation
  pinit (8) =  &array(0.,rtc.nWfs);
  pind (8) =  &amis;
  
  return [pinit,pind];
}

func mergesLayers(cn2h,alt,l0h,altsampling)
{

  step = (altsampling(dif))(1);
  cn2h2 = cn2h(sort(alt));
  alt2 = alt(sort(alt));
  l0h2 = l0h(sort(alt));
  //merges layers
  h2 = r2 = l2 = [];
  for(l=1;l<=numberof(altsampling);l++){
    if(is_array(alt2)){
      wi = where(alt2 >= altsampling(l)-step & alt2< altsampling(l)+step);
      if(is_array(wi)){
        hb = sum(alt2(wi)*cn2h2(wi))/sum(cn2h2(wi));
        h2 = grow(h2,hb);
        r2 = grow(r2,sum(cn2h2(wi)));
        l2 = grow(l2,min(l0h2(wi)));
        //remove the values that have been already taken into account
        wtmp = where((where4Vector(alt2,alt2(wi),notequal=1)));
        alt2 = alt2(wtmp);
        cn2h2 = cn2h2(wtmp);
        l0h2 = l0h2(wtmp)
      }
    }
  }
  perc = 100*cn2h/sum(cn2h);
  wh = where(alt > max(abs(altsampling)) & perc>=1);
  if(is_array(wh)){
    h2 = grow(h2,alt(wh));
    r2 = grow(r2,cn2h(wh));
    l2 = grow(l2,l0h(wh));
  }

  wdif = where(h2(dif)<dh);
  tmp2 = tmp = [r2,h2,l2];
    
    if(is_array(wdif)){
      for(d=1;d<=numberof(wdif);d++){
        cnh = tmp(wdif(d),1) + tmp(wdif(d)+1,1);
        alt = (tmp(wdif(d),2)*tmp(wdif(d),1) + tmp(wdif(d)+1,2)*tmp(wdif(d)+1,1) )/cnh;
        lh = tmp(wdif(d),3) + tmp(wdif(d)+1,3);
          
        tmp2(wdif(d),) = [cnh,alt,lh];
        tmp2(wdif(d)+1,) = 0;
      }
      tmp2 = tmp2(where(tmp2(,1,)!=0),);
    }
    
  return tmp2;
}

/*
 _                          
| |    ___  __ _ _ __ _ __  
| |   / _ \/ _` | '__| '_ \ 
| |__|  __/ (_| | |  | | | |
|_____\___|\__,_|_|  |_| |_|
                            
*/


func doLearn(Cdd,fitMethod,pinit,pind,Cnn,verb=,deriv=)
/* DOCUMENT fullMatrix = learn_lolo(fullCaa,fitMethod,pinit,pind,ptr_cn2h,Cnn=,verb=)

   Performs a learn fitting minimization on the matrix fullCaa using
   the initial guess given by pind. The parameters to be
   fitted is set thanks to pind. The algorithm uses the
   diagonal of the noise covariance matrix to be added to ftted
   covariance to regularize the noise.
 */

{
    
  if(is_void(fitMethod)){
    fitMethod = 1;
  }
  learn.transformationMatrix =  &identite(rtc.nSlopes);
  
  //TT management
  Cdd_ttr = handle_tilt_from_wfstype(Cdd,ttr = learn.ttr);
  Cnn_ttr = handle_tilt_from_wfstype(Cnn,ttr = learn.ttr);

  //Management of the diagonal of the matrix
  if(fitMethod == 2 ) {
    learn.diagonal = 0;
    // will subtract noise from diagonal of the matrix to fit only turbulence
    Cdd_ttr -= Cnn_ttr;
  }
  else {
    // case = "Do not fit diago" or "No noise at all"
    learn.diagonal = 1;
    takesDiag(Cdd_ttr) = 0.;
  }

  //transfers ptr init guess into tomo struct
  initStruct,pinit,pind;
  learn_fit.xshift   = 0;
  learn_fit.yshift   = 0;
  learn_fit.magnification = 0;
  learn_fit.theta    = 0;
  learn_fit.tracking = 0;
    
  fitEstim = packcoeffs( learn, indexFit );
  
  // .................. LEARN ....................//
  if( is_array(indexFit) ) {
    t = tic(10);
    if(fitMethod == 1){
      learn.runLearnTwoSteps = 0;
      res = lmfit_Learn( covMatModel, learn, fitEstim, Cdd_ttr, fit=indexFit, tol=1e-4,stdev=1,verb=verb,deriv=deriv);
    }else if(fitMethod == 2){
      if(learn.nl == 1){
        write,"\rI can't perform the 2 steps on a single layer turbulence, I do a classical Learn\r";
        learn.runLearnTwoSteps = 0;
        res = lmfit_Learn( covMatModel, learn, fitEstim, Cdd_ttr, fit=indexFit, tol=1e-4,stdev=1,verb=verb,deriv=deriv);
      }else{
        learn.runLearnTwoSteps = 1;
        res = Learn_2steps(learn,fitEstim,Cdd_ttr,indexFit,tol=1e-4,verb=verb,deriv=deriv);
      }
    }

    tinit = tac(10);
    if(verb)
      swrite(format="Done. Learn done in %f seconds",tinit);
    

    //...............TURBULENT PARAMETERS.....................//
   
    // updating tomo.learn structure 
    unpackcoeffs, fitEstim, cnh, alt,l0h,tracking,xshift,yshift,magnification,theta;
    ptrcn2h = [&abs(cnh), &alt, &abs(l0h),&tracking,&xshift,&yshift,&magnification,&theta,&tinit];
    if(is_array(res)){
      ptrcn2h = grow(ptrcn2h,[&res.chi2_last,&res.niter]);
      //getting uncertainties
      du = (*res.stdev)(indexFit);
      nl = learn.nl;
    }
    cnh2Struct, ptrcn2h;
    
    //Additionnal lm fit to get mis-registration

    if(learn.ttr == 1){
    //erase the learn_fit struct
    learn.ttr          = 0;
    learn_fit.cnh      = 0;
    learn_fit.altitude = 0;
    learn_fit.l0h      = 0;
    learn_fit.xshift   = *pind(5);
    learn_fit.yshift   = *pind(6);
    learn_fit.magnification   = *pind(7);
    learn_fit.theta    = *pind(8);
    learn_fit.tracking = [1,1,1];
    //forcing the diagonal to be null
    learn.diagonal = 1;
    takesDiag(Cdd) = 0; 
    fitEstim = packcoeffs(learn,indexFit);
    //tracking fitting
    res = lmfit_Learn( covMatModel, learn, fitEstim, Cdd, fit=indexFit, tol=1e-4,stdev=1,verb=verb,deriv=deriv);
    }
    
  //transfers results into learn struct
    unpackcoeffs, fitEstim, cnh, alt,l0h,tracking,xshift,yshift,magnification,theta;
    ptrcn2h = [&abs(cnh), &alt, &abs(l0h),&tracking,&xshift,&yshift,&magnification, &theta, &magnification,&tinit];
  cnh2Struct, ptrcn2h;
   if(is_array(res))
     ptrcn2h = grow(ptrcn2h,[&res.chi2_last,&res.niter]);
  }
  
  return ptrcn2h;  
}

func Learn_2steps(learn, &fitEstim, Cdd, indexFit,tol=,verb=,deriv=)
/* DOCUMENT  Learn_2steps(tomo, fitEstim, fullCaa, ww, indexFit)
   
   Performs the Learn 2 steps. First we subtract the ground layer from
   the covariance matrix of the measurements and the model. We fit the
   altitude and the strength of the altitude layers, depending on what
   we choose to fit before starting the Learn. Then, we fit the
   strength of the ground layer from the caa and the previous
   retrieved profile. In order to avoid some divergence of the fit, we
   have implemented two check of the divergence.
   
   The first one takes place after the first step and check if there
   is a layer with a strength greater than the previous computed R0
   (in meters) in the learn_do_it function. If it does, this diverging
   layer is fitted (altitude + strength) during the second step (from
   the entire covariance matrix) and the ground layer is let to 0.
   The second one takes place at the end of the second step and check
   if there is a layer greater than the ground layer. If it does, we
   fit the altitude and the strength of all layers whose strength are
   greater than the ground layer's, this one included.

   Olivier Martin
*/
{

  if(is_void(tol)) tol = 1e-4;
  
 
  //1. Computation of the transformation matrix for subtracting the ground layer
  P = ComputeTransMatrix2steps(2,ttr=learn.ttr);
  learn.transformationMatrix =  &P;//tomo.transformationMatrix used by calcMatCov

  //2. Subtract of the ground layer in the measurements
  Cdd_alt = (P(,+)*Cdd(,+))(,+)*P(+,);

  //3. Manage the tomo variables
  // If there's no layer with alt==0, then the lowest altitude is set to 0, and it
  // will NOT be fitted
  alt_array = learn.altitude(1:learn.nl);
  if(!sum(alt_array==0)){//if there is no ground layer in the initial guess

    where_altmin = where(alt_array == min(alt_array));//position of the lower altitude
    learn.altitude(where_altmin) = 0;//forces the minimal altitude to 0
    write, "\rThe minimal altitude was forced to 0 m\n",append=1;
  }

  indexGndLay = where(learn.altitude==0);

  //keep the desired fitting option in memory
  l0fitinit = learn_fit.l0h(indexGndLay);
  cnhfitinit = learn_fit.cnh(indexGndLay);
  
  learn.cnh(indexGndLay)= 0;//strength of the ground layer set to 0
  learn.tracking= 0;//tracking set to 0
  learn.xshift= learn.yshift = 0;
  learn_fit.cnh(indexGndLay) = 0;//no fit of the strength of the G.L
  learn_fit.altitude(indexGndLay) = 0;//no fit of the altitude of the G.L  
  learn_fit.l0h(indexGndLay) = 0;//no fit of L0 at the ground.
  learn_fit.tracking = 0;//no fit of the tracking error without G.L.
  learn_fit.xshift = learn_fit.yshift = 0;
  learn.diagonal = 1.;
  
  fitEstim = packcoeffs( learn, indexFit );
  //4. Fit without ground layer
  if(verb)
  write, " \rFirst step: fit of the altitude layers";
  takesDiag(Cdd_alt) = 0.;

  res = lmfit_Learn(covMatModel, learn, fitEstim, Cdd_alt, fit=indexFit, tol=tol,stdev=1,verb=verb,deriv=deriv );
  
  //5. Manage tomo variables
  unpackcoeffs, fitEstim, cnh, alt,l0h,tracking,xshift,yshift,magnification,theta;
  ptrcn2h = [&cnh, &alt, &l0h,&tracking,&xshift,&yshift,&magnification,&theta]; // will be saved later
  cnh2Struct, ptrcn2h;//update of the tomo struct

  learn_fit.altitude(1:learn.nl) = array(0.,learn.nl);//No altitude layer fitted
  learn_fit.cnh(1:learn.nl) = array(0.,learn.nl);//No strength fitted
  learn_fit.l0h(1:learn.nl) = array(0.,learn.nl);//No outerscale fitted
  
  
  learn.cnh(indexGndLay) = 20.0;
  learn.l0h(indexGndLay)= atm.L0;
  learn.tracking = [0,0,0];
  learn_fit.cnh(indexGndLay) = cnhfitinit;
  learn_fit.l0h(indexGndLay) = l0fitinit;
  
  fitEstim = packcoeffs( learn, indexFit );

  //6. Reinit of the transformation matrix
  learn.runLearnTwoSteps = 0;
  takesDiag(Cdd) = 0.;
  //7. Fit with GL
  if(verb)
    write, "\rSecond step: fit of the ground layers' strength\r";
  res = lmfit_Learn(covMatModel, learn, fitEstim, Cdd, fit=indexFit, tol=tol,stdev=1,verb=verb,deriv=deriv );
  
  return res;
}

/*
 ____  _                   _   
/ ___|| |_ _ __ _   _  ___| |_ 
\___ \| __| '__| | | |/ __| __|
 ___) | |_| |  | |_| | (__| |_ 
|____/ \__|_|   \__,_|\___|\__|
                               
                                                             _   
 _ __ ___   __ _ _ __   __ _  __ _  ___ _ __ ___   ___ _ __ | |_ 
| '_ ` _ \ / _` | '_ \ / _` |/ _` |/ _ \ '_ ` _ \ / _ \ '_ \| __|
| | | | | | (_| | | | | (_| | (_| |  __/ | | | | |  __/ | | | |_ 
|_| |_| |_|\__,_|_| |_|\__,_|\__, |\___|_| |_| |_|\___|_| |_|\__|
                             |___/                               
*/

func extendList(list, arr)
/* DOCUMENT
   concatenates a 1D array <arr> to an already existing 1D array
   <list>.
   <arr> is concatenated at the beginning of <list>, and the number of
   elements of <arr> is written at the end of <list>. Then the very
   last element list(0) is the number of concatenated arrays in
   <list>.
 */
{
  n = list(0);                 // number of elements stored in the list
  list = grow(arr, list(:-1));
  list = grow(list, numberof(arr), n+1);
  return list;
}


func grabList(&list)
/* DOCUMENT
   extracts a 1D array <arr> from an already existing 1D array
   <list>.
   <arr> is extracted at the beginning of <list>, and the number of
   elements of <arr> is written at list(-1). Then the very
   last element list(0) is the number of concatenated arrays in
   <list>.
   On output, <list> is modified so that the array <arr> is removed.
 */
{
  n = list(0);                 // number of elements stored in the list
  if( n>0 ) {
    p = long(list(-1));
    arr = list(1:p);
    list = list(p+1:-1);
    list(0) = n-1;
  } else {
    return [];
  }
  return arr;
}

func packcoeffs( learn , &indfit )
/* DOCUMENT

   Permet de grouper toutes les variables a fitter dans un seul
   tableau 1D. Ce tableau 1D sera "unpacked" par la fonction de calcul
   des mat de cov grace a la routine unpackcoeffs() qui extrait les
   differents tableaux en leur redonnant des noms plus explicites (r0,
   alt, etc..). Les fonctions packcoeffs() et unpackcoeffs() doivent
   absolument macher ensemble et etre symetriques l'une de l'autre.
   
   Pour rajouter un element a fitter:
     - le rajouter dans la structure tomo au bon endroit (variables.i)
     - penser a initialiser ces nouvelles variables, probablement dans
       initTomoParam() dans ce mm fichier plus haut
     - ajouter la ligne dans le 1er set de lignes de packcoeffs()
       (copier-coller 1 ligne ressemblante)
     - ajouter la ligne AU MM ENDROIT dans le 2ieme set (indfit_xxx ...)
     - ajouter le pendant, dans unpackcoeffs(), en prenant l'ordre INVERSE !!!
     - mettre bien a jout la liste des arguments de unpackcoeffs(), avec un
       nom explicite, et bien placé dans la liste
     - aller modifier l'appel a unpackcoeffs() dans les fonctions
         - calcMatCov() dans init_tomo_styc.i et modifier la fonctionnalite
       de la fonction
         - preview_targetpos() et learn_disp_coeffs() dans widget_learn.i
     - inserer les boites dans le glade avec des noms du genre variable_1_1
       et un handler=on_activate (et c'est tout y'a rien d'autre a modifier
       a ce niveau)
     - dans le learn, recuperer les contenus des nouvelles boites:
         - modifier la fn collectGuiVariables(self): dans widget_learn.py
         - modifier toggleFitBoolean() dans yorick/widget_learn.i
         - modifier updatePos_Gui2Tomo(..) dans yorick/widget_learn.i
         - modifier updatePos_Tomo2Gui(void) dans yorick/widget_learn.i
     - transferCn2h_to_tomoStruct( ptr ) doit aussi etre modifiee
     
   SEE ALSO:
 */
{
  nl = learn.nl;
  coeffList = [0.,0];
  coeffList = extendList(coeffList, learn.cnh(1:nl));
  coeffList = extendList(coeffList, learn.altitude(1:nl));
  coeffList = extendList(coeffList, learn.l0h(1:nl));
  coeffList = extendList(coeffList, learn.tracking);
  coeffList = extendList(coeffList, learn.xshift(1:rtc.nWfs));
  coeffList = extendList(coeffList, learn.yshift(1:rtc.nWfs));
  coeffList = extendList(coeffList, learn.magnification(1:rtc.nWfs));
  coeffList = extendList(coeffList, learn.theta(1:rtc.nWfs));
  
  indfit = [0.,0];
  indfit = extendList(indfit, learn_fit.cnh(1:nl));
  indfit = extendList(indfit, learn_fit.altitude(1:nl));
  indfit = extendList(indfit, learn_fit.l0h(1:nl));
  indfit = extendList(indfit, learn_fit.tracking);
  indfit = extendList(indfit, learn_fit.xshift(1:rtc.nWfs));
  indfit = extendList(indfit, learn_fit.yshift(1:rtc.nWfs));
  indfit = extendList(indfit, learn_fit.magnification(1:rtc.nWfs));
  indfit = extendList(indfit, learn_fit.theta(1:rtc.nWfs));

  indfit = indfit( :-long(indfit(0))-1 );   // cut the "tail" of the coeff pack
  indfit = where(indfit);
  
  return coeffList;
}


func unpackcoeffs(list,&cnh,&alt,&l0h,&tracking,&xshift,&yshift,&magnification,&theta)
/* DOCUMENT 
     Mode d'emploi : voir packcoeffs()
   SEE ALSO:
 */
{
  tmp = list;
  // they must be placed in REVERSE ORDER, because we UNstack here what has previously been stacked !!!
  theta         = grabList( tmp );
  magnification = grabList( tmp );
  yshift        = grabList( tmp );
  xshift        = grabList( tmp );
  tracking      = grabList( tmp );
  l0h           = grabList( tmp );
  alt           = grabList( tmp );
  cnh           = grabList( tmp );
}

func initStruct(pinit,pind)
/* DOCUMENT transfertInitGuessToTomoStruct,pinit,pind

   Transfers the initial guess and parameters to be fitted from
   pinit and pind into tomo struc
   
*/
{
  nl = learn.nl;
  //..............FITTED PARAMETERS........................//
  //cn2h
  (learn_fit.cnh)(1:nl) = *pind(1);
  //altitude
  (learn_fit.altitude)(1:nl) = *pind(2);
  //L0
  (learn_fit.l0h)(1:nl) = *pind(3);
  //tracking
  learn_fit.tracking = *pind(4);
  //pupil shift
  learn_fit.xshift(1:rtc.nWfs) = *pind(5);
  learn_fit.yshift(1:rtc.nWfs) = *pind(6);
  learn_fit.magnification(1:rtc.nWfs) = *pind(7);
  learn_fit.theta(1:rtc.nWfs) = *pind(8);

  //..............INITIAL GUESS........................//
  //cn2h
  (learn.cnh)(1:nl)      = *pinit(1);
  //altitude
  (learn.altitude)(1:nl) = *pinit(2);
  //L0
  (learn.l0h)(1:nl)      = *pinit(3);
  //tracking
  learn.tracking         = *pinit(4);
  //pupil shift
  learn.xshift(1:rtc.nWfs)        = *pinit(5);
  learn.yshift(1:rtc.nWfs)        = *pinit(6);
  learn.magnification(1:rtc.nWfs) = *pinit(7);
  learn.theta(1:rtc.nWfs)         = *pinit(8);
}

func cnh2Struct(ptrcnh)
/* DOCUMENT transfertInitGuessToTomoStruct,pinit,pind

   Transfers the initial guess and parameters to be fitted from
   pinit and pind into tomo struc
   
*/
{
  //Filling the learn struct
  nl = learn.nl;
  //cn2h
  (learn.cnh)(1:nl)      = *ptrcnh(1);
  //altitude
  (learn.altitude)(1:nl) = *ptrcnh(2);
  //L0
  (learn.l0h)(1:nl)      = *ptrcnh(3);
  //tracking
  learn.tracking         = *ptrcnh(4);
  //shift pupil;
  learn.xshift(1:rtc.nWfs)        = *ptrcnh(5);
  learn.yshift(1:rtc.nWfs)        = *ptrcnh(6);
  learn.magnification(1:rtc.nWfs) = *ptrcnh(7);
  learn.theta(1:rtc.nWfs)         = *ptrcnh(8);
}

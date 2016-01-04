/*
 _____ _ _   _   _             
|  ___(_) |_| |_(_)_ __   __ _ 
| |_  | | __| __| | '_ \ / _` |
|  _| | | |_| |_| | | | | (_| |
|_|   |_|\__|\__|_|_| |_|\__, |
                         |___/ 
*/
func fitCovarianceMatrix(dataarc,nl,ttr=,fitl0=,FitMethod=,fullHD=,tomores=,verb=,deriv=)
/* DOCUMENT ptrcn2h = fitCovarianceMatrix(*data.slopes_dis,3,fitl0=1,FitMethod=1,fullHD=0,tomores=5000,verb=1,deriv=0)

   Performs a quick fitting from input arguments.
 */
{
  data.learn.nl = nl;
  if(ttr) data.learn.ttr = 1;
  //.....INITIAL GUESS.....//
  tmp = initGuessLearn(nl,fitl0=fitl0,fullHD=fullHD,tomores=tomores);
  pinit = tmp(,1);
  pind = tmp(,2);

  //.....Compute covariance matrix of turbulence with square ssp.....//
  dataarc -= dataarc(,avg);
  Cdd = dataarc(,+)*dataarc(,+)/dimsof(dataarc)(0);

  //...... Derivating the noise covariance matrix ......//
  Cnn = 0*Cdd;
  takesDiag(Cnn)  = getNoiseVar(dataarc);

  //fit
  if(is_void(FitMethod)){
    FitMethod = 1;
  }

  ptrcn2h = learn(Cdd,FitMethod,pinit,pind,Cnn,verb=verb,deriv=deriv);
 
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
  D = data.tel.diam;
  //Base between each stars
  for(i=1;i<=n-1;i++){
    if(useTS || i!=n+1){
      Btmp = sqrt((posX(i) - posX(i+1:))^2 + (posY(i) - posY(i+1:))^2);
      if(alt_lgs(i)){
        Blim = D*206265/alt_lgs(i);
        Btmp*=(Btmp > Blim);//if altitude max > LGS altitude: but baseline to 0
      }
      B = grow(B,Btmp);
    }
  }
  B = B(where(B!=0))/206265.;
  
  if(anyof(B)){
    Bmax = max(B);Bmin=min(B);
    if(altituderes)
      return  D/Bmax/data.wfs(data.its).sX;
    else 
      return D/Bmin;
  }

}

func initGuessLearn(nbl,fitl0=,fullHD=,tomores=)
  /* DOCUMENT [ptr_initarray,pind] = initGuessLearn(nb_layers,posX,posY,fullHD=,tomores=)

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
  altmax = findAltitudeMax(data.wfs.x,data.wfs.y,useTS=1);//find the maximal altitude seen by Canary without TS
  pinit = pind = array(pointer,6);

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
  l0h = data.turbu.L0;
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
  if(data.learn.ttr == 0){
    pind(4) = &array(1,3);
  }else{pind(4) = &array(0,3);}

  //shift pupil
  pinit (5) = pinit (6) = &array(0.,data.nwfs);
  pind (5) = pind (6) = &array(0.,data.nwfs);
  
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


func learn(Cdd,fitMethod,pinit,pind,Cnn,verb=,deriv=)
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
  data.learn.transformationMatrix =  &identite(data.Nslopes);

  //TT management
  Cdd_ttr = handle_tilt_from_wfstype(Cdd,ttr = data.learn.ttr);
  Cnn_ttr = handle_tilt_from_wfstype(Cnn,ttr = data.learn.ttr);

  //Management of the diagonal of the matrix
  if(fitMethod == 2 ) {
    data.learn.diagonal = 0;
    // will subtract noise from diagonal of the matrix to fit only turbulence
    Cdd_ttr -= Cnn_ttr;
  }
  else {
    // case = "Do not fit diago" or "No noise at all"
    data.learn.diagonal = 1;
    takesDiag(Cdd_ttr) = 0.;
  }

  //transfers ptr init guess into tomo struct
  initStruct,pinit,pind;
  fitEstim = packcoeffs( data, indexFit );
  
  // .................. LEARN ....................//
  if( is_array(indexFit) ) {
    t = tic(10);
    if(fitMethod == 1){
      data.learn.runLearnTwoSteps = 0;
      res = lmfit_Learn( covMatModel, data, fitEstim, Cdd_ttr, fit=indexFit, tol=1e-4,stdev=1,verb=verb,deriv=deriv);
    }else if(fitMethod == 2){
      if(data.learn.nl == 1){
        write,"\rI can't perform the 2 steps on a single layer turbulence, I do a classical Learn\r";
        data.learn.runLearnTwoSteps = 0;
        res = lmfit_Learn( covMatModel, data, fitEstim, Cdd_ttr, fit=indexFit, tol=1e-4,stdev=1,verb=verb,deriv=deriv);
      }else{
        data.learn.runLearnTwoSteps = 1;
        res = Learn_2steps(data,fitEstim,Cdd_ttr,indexFit,tol=1e-4,verb=verb,deriv=deriv);
      }
    }

    tinit = tac(10);
    if(verb)
      swrite(format="Done. Learn done in %f seconds",tinit);
    

    //...............TURBULENT PARAMETERS.....................//
   
    // updating tomo.learn structure 
    unpackcoeffs, fitEstim, cnh, alt,l0h,tracking,xshift,yshift;
    ptrcn2h = [&abs(cnh), &alt, &abs(l0h),&tracking,&xshift,&yshift,&tinit];
    if(is_array(res)){
      ptrcn2h = grow(ptrcn2h,[&res.chi2_last,&res.niter]);
      //getting uncertainties
      du = (*res.stdev)(indexFit);
      nl = data.learn.nl;
      data.uncertainties.l0h(1:nl) = du(1:nl);
      data.uncertainties.altitude(1:nl) = du(nl+1:2*nl);
      data.uncertainties.cnh(1:nl) = du(2*nl+1:3*nl);
    }
    cnh2Struct, ptrcn2h;
    
  //..... Additionnal lm fit to get the tracking if it was not already done .....//
  if(data.learn.ttr == 1){
    //erase the learn_fit struct
    data.learn.ttr = 0;
    data.learn_fit.cnh = 0;
    data.learn_fit.altitude = 0;
    data.learn_fit.l0h=0;
    data.learn_fit.xshift = 0;
    data.learn_fit.yshift = 0;
    data.learn_fit.tracking = [1,1,1];
    //forcing the diagonal to be null
    data.learn.diagonal = 1;
    takesDiag(Cdd) = 0; 
    fitEstim = packcoeffs(data,indexFit);
    //tracking fitting
    res = lmfit_Learn( covMatModel, data, fitEstim, Cdd, fit=indexFit, tol=1e-4,stdev=1,verb=verb,deriv=deriv);
    data.uncertainties.tracking = (*res.stdev)(where(*res.stdev !=0));
  }
  
  //transfers results into learn struct
  unpackcoeffs, fitEstim, cnh, alt,l0h,tracking,xshift,yshift;
  ptrcn2h = [&abs(cnh), &alt, &abs(l0h),&tracking,&xshift,&yshift,&tinit];
  cnh2Struct, ptrcn2h;
  }
  
  return ptrcn2h;  
}

func Learn_2steps(data, &fitEstim, Cdd, indexFit,tol=,verb=,deriv=)
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
  P = ComputeTransMatrix2steps(2,ttr=data.learn.ttr);
  data.learn.transformationMatrix =  &P;//tomo.transformationMatrix used by calcMatCov

  //2. Subtract of the ground layer in the measurements
  Cdd_alt = (P(,+)*Cdd(,+))(,+)*P(+,);

  //3. Manage the tomo variables
  // If there's no layer with alt==0, then the lowest altitude is set to 0, and it
  // will NOT be fitted
  alt_array = data.learn.altitude(1:data.learn.nl);
  if(!sum(alt_array==0)){//if there is no ground layer in the initial guess

    where_altmin = where(alt_array == min(alt_array));//position of the lower altitude
    data.learn.altitude(where_altmin) = 0;//forces the minimal altitude to 0
    write, "\rThe minimal altitude was forced to 0 m\n",append=1;
  }

  indexGndLay = where(data.learn.altitude==0);

  //keep the desired fitting option in memory
  l0fitinit = data.learn_fit.l0h(indexGndLay);
  cnhfitinit = data.learn_fit.cnh(indexGndLay);
  
  data.learn.cnh(indexGndLay)= 0;//strength of the ground layer set to 0
  data.learn.tracking= 0;//tracking set to 0
  data.learn.xshift= data.learn.yshift = 0;
  data.learn_fit.cnh(indexGndLay) = 0;//no fit of the strength of the G.L
  data.learn_fit.altitude(indexGndLay) = 0;//no fit of the altitude of the G.L  
  data.learn_fit.l0h(indexGndLay) = 0;//no fit of L0 at the ground.
  data.learn_fit.tracking = 0;//no fit of the tracking error without G.L.
  data.learn_fit.xshift = data.learn_fit.yshift = 0;
  data.learn.diagonal = 1.;
  
  fitEstim = packcoeffs( data, indexFit );
  //4. Fit without ground layer
  write, "\rFirst step: fit of the altitude layers";
  takesDiag(Cdd_alt) = 0.;

  res = lmfit_Learn(covMatModel, data, fitEstim, Cdd_alt, fit=indexFit, tol=tol,stdev=1,verb=verb,deriv=deriv );
  
  //5. Manage tomo variables
  unpackcoeffs, fitEstim, cnh, alt,l0h,tracking,xshift,yshift;
  ptrcn2h = [&cnh, &alt, &l0h,&tracking,&xshift,&yshift]; // will be saved later
  cnh2Struct, ptrcn2h;//update of the tomo struct

  data.learn_fit.altitude(1:data.learn.nl) = array(0.,data.learn.nl);//No altitude layer fitted
  data.learn_fit.cnh(1:data.learn.nl) = array(0.,data.learn.nl);//No strength fitted
  data.learn_fit.l0h(1:data.learn.nl) = array(0.,data.learn.nl);//No outerscale fitted
  
  
  data.learn.cnh(indexGndLay) = 20.0;
  data.learn.l0h(indexGndLay)= data.turbu.L0;
  data.learn.tracking = [0,0,0];
  data.learn_fit.cnh(indexGndLay) = cnhfitinit;
  data.learn_fit.l0h(indexGndLay) = l0fitinit;
  
  fitEstim = packcoeffs( data, indexFit );

  //6. Reinit of the transformation matrix
  data.learn.runLearnTwoSteps = 0;
  takesDiag(Cdd) = 0.;
  //7. Fit with GL
  write, "\rSecond step: fit of the ground layers' strength\r";
  res = lmfit_Learn(covMatModel, data, fitEstim, Cdd, fit=indexFit, tol=tol,stdev=1,verb=verb,deriv=deriv );
  
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

func packcoeffs( data , &indfit )
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
  nl = data.learn.nl;
  coeffList = [0.,0];
  coeffList = extendList(coeffList, data.learn.cnh(1:nl));
  coeffList = extendList(coeffList, data.learn.altitude(1:nl));
  coeffList = extendList(coeffList, data.learn.l0h(1:nl));
  coeffList = extendList(coeffList, data.learn.tracking);
  coeffList = extendList(coeffList, data.learn.xshift);
  coeffList = extendList(coeffList, data.learn.yshift);
  
  indfit = [0.,0];
  indfit = extendList(indfit, data.learn_fit.cnh(1:nl));
  indfit = extendList(indfit, data.learn_fit.altitude(1:nl));
  indfit = extendList(indfit, data.learn_fit.l0h(1:nl));
  indfit = extendList(indfit, data.learn_fit.tracking);
  indfit = extendList(indfit, data.learn_fit.xshift);
  indfit = extendList(indfit, data.learn_fit.yshift);

  indfit = indfit( :-long(indfit(0))-1 );   // cut the "tail" of the coeff pack
  indfit = where(indfit);
  
  return coeffList;
}


func unpackcoeffs(list,&cnh, &alt,&l0h,&tracking,&xshift,&yshift)
/* DOCUMENT 
     Mode d'emploi : voir packcoeffs()
   SEE ALSO:
 */
{
  tmp = list;
  // WAAAAAAAAARNING !!
  //
  // they must be placed in REVERSE ORDER, because we UNstack here what has previously been stacked !!!
  xshift = grabList( tmp);
  yshift = grabList (tmp);
  tracking = grabList( tmp );
  l0h = grabList( tmp );
  alt = grabList( tmp );
  cnh = grabList( tmp );
}

func initStruct(pinit,pind)
/* DOCUMENT transfertInitGuessToTomoStruct,pinit,pind

   Transfers the initial guess and parameters to be fitted from
   pinit and pind into tomo struc
   
*/
{
  nl = data.learn.nl;
  //..............FITTED PARAMETERS........................//
  //cn2h
  (data.learn_fit.cnh)(1:nl) = *pind(1);
  //altitude
  (data.learn_fit.altitude)(1:nl) = *pind(2);
  //L0
  (data.learn_fit.l0h)(1:nl) = *pind(3);
  //tracking
  data.learn_fit.tracking = *pind(4);
  //pupil shift
  data.learn_fit.xshift = *pind(5);
  data.learn_fit.yshift = *pind(6);

  //..............INITIAL GUESS........................//
  //cn2h
  (data.learn.cnh)(1:nl) = *pinit(1);
  //altitude
  (data.learn.altitude)(1:nl) = *pinit(2);
  //L0
  (data.learn.l0h)(1:nl) = *pinit(3);
  //tracking
  data.learn.tracking = *pinit(4);
  //pupil shift
  data.learn.xshift = *pinit(5);
  data.learn.yshift = *pinit(6);
}

func cnh2Struct(ptrcnh)
/* DOCUMENT transfertInitGuessToTomoStruct,pinit,pind

   Transfers the initial guess and parameters to be fitted from
   pinit and pind into tomo struc
   
*/
{
  nl = data.learn.nl;
  //cn2h
  (data.learn.cnh)(1:nl) = *ptrcnh(1);
  //altitude
  (data.learn.altitude)(1:nl) = *ptrcnh(2);
  //L0
  (data.learn.l0h)(1:nl) = *ptrcnh(3);
  //tracking
  data.learn.tracking = *ptrcnh(4);
  //shift pupil;
  data.learn.xshift = *ptrcnh(5);
  data.learn.yshift = *ptrcnh(6);
}

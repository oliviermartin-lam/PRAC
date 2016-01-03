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
  
  indfit = [0.,0];
  indfit = extendList(indfit, data.learn.cnhfit(1:nl));
  indfit = extendList(indfit, data.learn.altfit(1:nl));
  indfit = extendList(indfit, data.learn.l0hfit(1:nl));

  indfit = indfit( :-long(indfit(0))-1 );   // cut the "tail" of the coeff pack
  indfit = where(indfit);
  
  return coeffList;
}


func unpackcoeffs(list,&cnh, &alt,&l0h)
/* DOCUMENT 
     Mode d'emploi : voir packcoeffs()
   SEE ALSO:
 */
{
  tmp = list;
  // WAAAAAAAAARNING !!
  //
  // they must be placed in REVERSE ORDER, because we UNstack here what has previously been stacked !!!
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
  (data.learn.cnhfit)(1:nl) = *pind(1);
  //altitude
  (data.learn.altfit)(1:nl) = *pind(2);
  //L0
  (data.learn.l0hfit)(1:nl) = *pind(3);

  //..............INITIAL GUESS........................//
  //cn2h
  (data.learn.cnh)(1:nl) = *pinit(1);
  //altitude
  (data.learn.altitude)(1:nl) = *pinit(2);
  //L0
  (data.learn.l0h)(1:nl) = *pinit(3);
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
}

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
  
  //Management of the diagonal of the matrix
  if(fitMethod == 2 ) {
    data.learn.diagonal = 0;
    // will subtract noise from diagonal of the matrix to fit only turbulence
    Cdd -= Cnn;
  }
  else {
    // case = "Do not fit diago" or "No noise at all"
    data.learn.diagonal = 1;
    takesDiag(Cdd) = 0.;
  }

  //transfers ptr init guess into tomo struct
  initStruct,pinit,pind;
  fitEstim = packcoeffs( data, indexFit );
  // .................. LEARN ....................//
  if( is_array(indexFit) ) {
    t = tic(10);
    if(fitMethod == 1){
      data.learn.runLearnTwoSteps = 0;
      res = lmfit_Learn( covMatModel, data, fitEstim, Cdd, fit=indexFit, tol=1e-4,stdev=1,verb=verb,deriv=deriv);
    }else if(fitMethod == 2){
      if(data.learn.nl == 1){
        write,"\rI can't perform the 2 steps on a single layer turbulence, I do a classical Learn\r";
        data.learn.runLearnTwoSteps = 0;
        res = lmfit_Learn( covMatModel, data, fitEstim, Cdd, fit=indexFit, tol=1e-4,stdev=1,verb=verb,deriv=deriv);
      }else{
        data.learn.runLearnTwoSteps = 1;
        res = Learn_2steps(data,fitEstim,Cdd,indexFit,tol=1e-4,verb=verb,deriv=deriv);
      }
    }

    tinit = tac(10);
    swrite(format="Done. Learn done in %f seconds",tinit);
    
    //...............TURBULENT PARAMETERS.....................//
   
    // updating tomo.learn structure 
    unpackcoeffs, fitEstim, cnh, alt,l0h;
    ptrcn2h = [&abs(cnh), &alt, &abs(l0h),&tinit];
    if(is_array(res)){
      ptrcn2h = grow(ptrcn2h,[&res.chi2_last,&res.niter]);
    }
    //transfers results into tomo struc
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
  P = ComputeTransMatrix2steps(2);
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
  l0fitinit = data.learn.l0hfit(indexGndLay);
  cnhfitinit = data.learn.cnhfit(indexGndLay);
  
  data.learn.cnh(indexGndLay)= 0;//strength of the ground layer set to 0
  data.learn.l0h(indexGndLay)= 0;//outerscale of the ground layer set to 0
  data.learn.cnhfit(indexGndLay) = 0;//no fit of the strength of the G.L
  data.learn.altfit(indexGndLay) = 0;//no fit of the altitude of the G.L  
  data.learn.l0hfit(indexGndLay) = 0;//no fit of L0 without G.L.
  
  data.learn.diagonal = 1.;
  fitEstim = packcoeffs( data, indexFit );
  //4. Fit without ground layer
  write, "\rFirst step: fit of the altitude layers";
  takesDiag(Cdd_alt) = 0.;

  res = lmfit_Learn(covMatModel, data, fitEstim, Cdd_alt, fit=indexFit, tol=tol,stdev=1,verb=verb,deriv=deriv );
  
  //5. Manage tomo variables
  unpackcoeffs, fitEstim, cnh, alt,l0h;
  ptrcn2h = [&cnh, &alt, &l0h]; // will be saved later
  cnh2Struct, ptrcn2h;//update of the tomo struct

  data.learn.altfit(1:data.learn.nl) = array(0.,data.learn.nl);//No altitude layer fitted
  data.learn.cnhfit(1:data.learn.nl) = array(0.,data.learn.nl);//No strength fitted
  data.learn.l0hfit(1:data.learn.nl) = array(0.,data.learn.nl);//No outerscale fitted
  
  
  data.learn.cnh(indexGndLay) = 20.0;
  data.learn.l0h(indexGndLay)= data.turbu.L0;
  data.learn.cnhfit(indexGndLay) = cnhfitinit;
  data.learn.l0hfit(indexGndLay) = l0fitinit;
  
  fitEstim = packcoeffs( data, indexFit );

  //6. Reinit of the transformation matrix
  data.learn.runLearnTwoSteps = 0;
  takesDiag(Cdd) = 0.;
  //7. Fit with GL
  write, "\rSecond step: fit of the ground layers' strength\r";
  res = lmfit_Learn(covMatModel, data, fitEstim, Cdd, fit=indexFit, tol=tol,stdev=1,verb=verb,deriv=deriv );
  
  return res;
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
  pinit = pind = array(pointer,3);

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
    l0h = grow(data.turbu.L0,array(100.,nbl-1));
  }else if(fitl0 !=2){
    l0h = grow(data.turbu.L0,array(1e5,nbl-1));
  }
  pinit(3) = &l0h;
  if(is_void(fitl0)| fitl0 == 0){
    pind(3) = &(array(0.,nbl));
  }else if(fitl0 == 1){
    pind(3) = &(grow(1,array(0.,nbl-1)));
  }else if(fitl0 == 2){
    pind(3) = &(array(1.,nbl));
  }
  

  
  return [pinit,pind];
}

func fitCovarianceMatrix(dataarc,nl,fitl0=,FitMethod=,fullHD=,tomores=,verb=,deriv=)
/* DOCUMENT ptrcn2h = fitCovarianceMatrix(*data.slopes_dis,3,fitl0=1,FitMethod=1,fullHD=0,tomores=5000,verb=1,deriv=0)

   Performs a quick fitting from input arguments.
 */
{
  data.learn.nl = nl;
  //.....INITIAL GUESS.....//
  tmp = initGuessLearn(nl,fitl0=fitl0,fullHD=fullHD,tomores=tomores);
  pinit = tmp(,1);
  pind = tmp(,2);

  //.....Compute covariance matrix of turbulence with square ssp.....//
  dataarc -= dataarc(,avg);
  Cdd = dataarc(,+)*dataarc(,+)/dimsof(dataarc)(0);
  Cdd = handle_tilt_from_wfstype( Cdd );

  //...... Derivating the noise covariance matrix ......//
  Cnn = 0*Cdd;
  takesDiag(Cnn)  = noise_inPix2(dataarc);
  Cnn = handle_tilt_from_wfstype( Cnn );

  //fit
  if(is_void(FitMethod)){
    FitMethod = 1;
  }
  ptrcn2h = learn(Cdd,FitMethod,pinit,pind,Cnn,verb=verb,deriv=deriv);
 
  return ptrcn2h;
}







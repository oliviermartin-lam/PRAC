include, "pracMain.i";

func testForAll(fct,verb=)
/* DOCUMENT r = testForAll(testFittingOtf,verb=1);

 */
{
  include, "pracConfig.i",1;
  //finding directories
  listdir = (listFile(dataDir))(:-1);
  ndir = numberof(listdir);
  res = [];

  //managing directories
  for(j=1;j<=ndir;j++){
    dir      = listdir(j);
    pathdata =  dir + "/";
    listtl =  listFile(dataDir + dir + "/slopestl/");
    w = where(strfind("script",listtl)(0,) != -1);
    listtl = listtl(w);
    ntl    = numberof(listtl);
    procDir = dir + "/";
    write,format="Processing directory %s\n",dir;

    //loop on file
    for(i=1;i<=ntl;i++){
      //determining the airmass
      timedata = strpart(listtl(i),21:29);
      res = grow(res,fct(timedata,Dir=procDir,verb=0));
      write,format="Job done: %.3g%s\r",100.*i/ntl,"%";
    }
  }
  return res;
}

func testFittingOtf(timedata,Dir=,verb=)
/* DOCUMENT testFittingOtf("02h30m48s");

 */
{
  //Instantiate the required structures
  include, "pracConfig.i",1;
  if(Dir){
    procDir = Dir;
    dataDirRoot = dataDir + procDir;
  }

  if(verb)
    write,format = "%s\r","Initializing the required structures";
  include, "pracStructConfig.i",1;
  define_structs,timedata;


  r0tot = atm.r0 *(cam.lambda/atm.lambda)^(1.2);
  
  // Teescope OTF
  if(verb)
    write,format = "%s\r","Computing the telescope OTF";
  
  OTF_tel = OTF_telescope(tel.diam,tel.obs,tel.nPixels,tel.pixSize);
   
  //Determine the fitting OTF using three different frequency mask

  if(verb)
    write,format = "%s\r","Computing the square fitting OTF";
  otfSquare = computeOTFfitting("square");
  if(verb)
    write,format = "%s\r","Computing the circle fitting OTF";
  otfCircle = computeOTFfitting("circle");
  if(verb)
    write,format = "%s\r","Computing the influence fitting OTF";
  otfInflue = computeOTFfitting("influence");

  //Compute the fitting error from calibrated relation E2E-based
  sigFit = computesFittingError(tel.diam,r0tot,atm.lambda*1e9);

  //Determine Strehl ratios

  srSquare = sum(otfSquare * OTF_tel);
  srCircle = sum(otfCircle * OTF_tel);
  srInflue = sum(otfInflue * OTF_tel);
  srE2E    = exp(-(sigFit*1e-9*2*pi/cam.lambda)^2);

  return [srSquare,srCircle,srInflue,srE2E];
}

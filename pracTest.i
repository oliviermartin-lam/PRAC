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


func plotTFMOAO(void){

  N = 100;
  tret = 3e-3;
  Fe = 150.;
  BP = 500;
  f = span(0,Fe,N);
  g = [0.1,0.5,1.];
  ng = numberof(g);
 

  winkill,0;window,0,dpi=90,style="aanda.gs";
  logxy,1,0;
  for(i=1;i<=numberof(g);i++){
    hcor = hcorMoao(f,Fe,tret,g(i),BP);
    hcor = 20*log10(abs(hcor(2:N/2)));
    plg,hcor,f(2:N/2),marks=0;
    x = Fe/10.; y = hcor(int(N/10.));
    plt,"gain = " + var2str(g(i)),x,y,tosys=1;
  }
  plg,0*f(2:N/2),f(2:N/2),type=2,marks=0;
  range,-30,10;
  xytitles,"Frequency (Hz)","Amplitude (dB)";
  pltitle, "MOAO correction transfer function";
  
  //TF bruit
  
  winkill,1;window,1,dpi=90,style="aanda.gs";
  logxy,1,0;
  for(i=1;i<=numberof(g);i++){
    hbo = hboMoao(f,Fe,tret,g(i),BP);
    hbo = 20*log10(abs(hbo(2:N/2)));
    plg,hbo,f(2:N/2),marks=0;
    x = Fe/10.; y = hbo(int(N/10.));
    plt,"gain = " + var2str(g(i)),x,y,tosys=1;
  }
  
  plg,0*f(2:N/2),f(2:N/2),type=2,marks=0;
  range,-30,10;
  xytitles,"Frequency (Hz)","Amplitude (dB)";
  pltitle, "MOAO noise rejection transfer function";
}

func hcreep(freq,c,Fe)
{
  Te = 1./Fe;
  p = 2i*pi*freq + 1e-12;
  z = exp(p*Te);
  euler = -.0577;
  
  hcreep = (z*z - c*euler*(z-1)-c*(z-1)*log(z))/(z*(z-1));

  return hcreep;
}

func creepTransferFunction(void)
/* DOCUMENT

 */

{

N = 100;
  tret = 3e-3;
  Fe = 150.;
  BP = 500;
  f = span(0,Fe,N);
  g = [0.1,0.5,1.];
  ng = numberof(g);
 

  winkill,0;window,0,dpi=90,style="aanda.gs";
  logxy,1,0;
  for(i=1;i<=numberof(g);i++){
    hcor = hcreep(f,c);
    hcor = 20*log10(abs(hcor(2:N/2)));
    plg,hcor,f(2:N/2),marks=0;
    x = Fe/10.; y = hcor(int(N/10.));
    plt,"gain = " + var2str(g(i)),x,y,tosys=1;
  }
  plg,0*f(2:N/2),f(2:N/2),type=2,marks=0;
  range,-30,10;
  xytitles,"Frequency (Hz)","Amplitude (dB)";
  pltitle, "Creeping transfer function";

}

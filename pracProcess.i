include, "pracMain.i";




func plotProfiles(void)
/* DOCUMENT plotProfiles;
  
 */
{

  cnh = readfits("results/seeingProfile.fits");
  l0h = readfits("results/l0hProfile.fits");
  vh = readfits("results/vhProfile.fits");
  alt = readfits("results/altProfile.fits");
  cnh0 = readfits("results/seeing.fits");
  cnhg = readfits("results/groundSeeing.fits");
  cnhalt = readfits("results/altSeeing.fits");
  l0eff = readfits("results/effL0.fits");
  l0g = readfits("results/groundL0.fits");
  l0alt = readfits("results/altL0.fits");
  v0 = readfits("results/effV0.fits");
  v0g = readfits("results/groundV0.fits");
  v0alt = readfits("results/altV0.fits");
  
  if(0){
    //initialiazing vectors...
    cnh = alt = l0h = vh = dirh = [];//fitted values
    dcnh = dalt = dl0h = dvh = ddirh = [];//uncertainties

    listdir = listFile("profiles/");
    ndir = numberof(listdir);
    //managing directories
    for(j=1;j<=ndir;j++){
      dir      = listdir(j);
      pathdata = "profiles/" + dir + "/";
      listtl = listFile(pathdata);
      ntl    = numberof(listtl);
      dataDirRoot = dataDir + dir + "/";
      write,format="Processing directory %s\n",dir;
      //loop on file
      for(i=1;i<=ntl;i++){
        //determining the airmass
        timetl = strpart(listtl(i),21:29);
        restorefits,"slopestl",timetl,pathtl,fake=1;
        object = strcase(1,readFitsKey(pathtl,"OBJECT"));
        tmp = givesCoordinatesFromObject(object);
        array_ra = str2flt(decoupe(tmp(1),' '));
        array_dec = str2flt(decoupe(tmp(2),' '));
        airm = airmassFromDate(extractDate(pathtl),array_ra,array_dec);
        //growing results
        ptrprof = readfits(pathdata + listtl(i));
        cnh     = grow(cnh,*ptrprof(1)/airm);
        alt     = grow(alt,*ptrprof(2)/airm);
        l0h     = grow(l0h,*ptrprof(3));
        track   = grow(track,*ptrprof(4));
        vh      = grow(vh,*ptrprof(5));
        dirh    = grow(dirh,*ptrprof(6));
        dcnh    = grow(dcnh,*ptrprof(7)/airm);
        dalt    = grow(dalt,*ptrprof(8)/airm);
        dl0h    = grow(dl0h,*ptrprof(9));
        dtrack  = grow(dtrack,*ptrprof(10));
        dvh     = grow(dvh,*ptrprof(11));
        ddirh   = grow(ddirh,*ptrprof(12));
      
        write,format="Loading file:%.3g%s\r",100.*i/ntl,"%";
      
      }
    }
  
    
    alt    /= 1000;
    //computing integrated value
    cnh0 = l0eff = v0 = cnhalt = l0alt = v0alt = cnhg = l0g = v0g = [];
    for(i=1;i<=numberof(cnh)/5;i++){
      cnhi   = cnh(1+(i-1)*5:5*i);
      cnha   = cnhi(2:);
      l0hi   = l0h(1+(i-1)*5:5*i);
      wi     = where(l0hi<=100.)
        vhi    = vh(1+(i-1)*5:5*i);
      //global values
      cnh0   = grow(cnh0,0.103/sum(cnhi)^(-3/5.));
      if(is_array(wi))
        l0eff  = grow(l0eff,(sum(l0hi(wi)^(5/3.)*cnhi(wi))/sum(cnhi(wi)))^(3/5.));
      v0     = grow(v0,(sum(cnhi*vhi^(5/3.))/sum(cnhi))^(3/5.));
      //ground values
      cnhalt = grow(cnhalt,0.103/sum(cnha)^(-3/5.));
      if(numberof(wi) == 1 && wi(1) != 1){
        l0alt  = grow(l0alt,(sum(cnhi(wi)*l0hi(wi)^(5/3.))/sum(cnhi(wi)))^(3/5.));
      }else if(is_array(wi) && numberof(wi)!=1){
        l0alt  = grow(l0alt,(sum(cnhi(wi(2:))*l0hi(wi(2:))^(5/3.))/sum(cnhi(wi(2:))))^(3/5.));
      }
      v0alt  = grow(v0alt,(sum(cnha*vhi(2:)^(5/3.))/sum(cnha))^(3/5.));
      //values integrated onto altitude
      cnhg   = grow(cnhg,0.103/cnhi(1)^(-3/5.));
      l0g    = grow(l0g,l0hi(1));
      v0g    = grow(v0g,vhi(1));
    }

    cnh     = 0.103/cnh^(-3/5.);
    dcnh    = 0.103*3/5.*dcnh/cnh^(2/5.); 
    wa      = where(alt>=20);
    alt(wa) = 20.;
    w0 = where(l0h<=40);
    write,format="%.3g%s of jerk outer scale\n",100.*(1.- double(numberof(w0))/numberof(l0h)),"%s";

    writefits,"seeingProfile.fits",cnh;
    writefits,"l0hProfile.fits",l0h;
    writefits,"vhProfiles.fits",vh;
    writefits,"altProfile.fits",alt;
    writefits,"seeing.fits",cnh0;
    writefits,"groundSeeing.fits",cnhg;
    writefits,"altSeeing.fits",cnhalt;
    writefits,"effL0.fits",l0eff;
    writefits,"groundL0.fits",l0g;
    writefits,"altL0.fits",l0alt;
    writefits,"effV0.fits",v0;
    writefits,"groundV0.fits",v0g;
    writefits,"altV0.fits",v0alt;
  }
  write,format="%d processed files\n",numberof(cnh)/5;
  //..... DISPLAY .....//
  
  //Plotting altitude versus cnh
  winkill,0;window,0,dpi=90,style="aanda.gs";
  plmk,alt,cnh,marker=4,msize=0.1;
  range,-1,25;
  limits,0,2.;
  xytitles,"Seeing at 500 nm (arcsec)","Altitude (km)";
  pdf, "cnh_profile.pdf";

  //Plotting seeing histograms
  winkill,4;window,4,dpi=90,style="aanda.gs";
  res = .04;cnhlim = 2.;
  ns = int(minmax(cnh0(where(cnh0<cnhlim)))(dif)(1)/res-1.);
  histo,cnh0(where(cnh0<cnhlim)),ns;
  histo,cnhg(where(cnhg<cnhlim)),ns;
  histo,cnhalt(where(cnhalt<cnhlim)),int(ns/1.05);
  xytitles,"Seeing at 500 nm (arcsec)","Counts";
  lim = limits();
  plt,"A: Global seeing",lim(2)/2,lim(4)/1.5,tosys=1;
  plt,"B: Ground value ",lim(2)/2,lim(4)/1.5-lim(4)/8,tosys=1;
  plt,"C: Altitude value",lim(2)/2,lim(4)/1.5-lim(4)/4,tosys=1;
  limits,-.1,2.2;
  pdf, "cnh_histo.pdf";
  w0 = where(cnh0<=cnhlim);
  write,format="Average/1-sigma global seeing         = %.3g''/ %.3g''\n",avg(cnh0(w0)),(cnh0(w0))(rms);
  w0 = where(cnhg<=cnhlim);
  write,format="Average/1-sigma ground seeing         = %.3g''/ %.3g''\n",avg(cnhg(w0)),(cnhg(w0))(rms);
  w0 = where(cnhalt<=cnhlim);
  write,format="Average/1-sigma altitude seeing       = %.3g''/ %.3g''\n",avg(cnhalt(w0)),(cnhalt(w0))(rms);
  
  //Plotting altitude versus l0h
  winkill,1;window,1,dpi=90,style="aanda.gs";
  w0 = where(l0h<=100.);
  plmk,alt(w0),l0h(w0),marker=4,msize=0.1;
  range,-1,25;
  limits,-5,100;
  xytitles,"L_0_(h) (m)","Altitude (km)";
  pdf, "L0h_profile.pdf";

  //Plotting outer scale histograms
  winkill,5;window,5,dpi=90,style="aanda.gs";
  res = 1.;l0lim=100.;
  ns = int(minmax(l0eff(where(l0eff<=l0lim)))(dif)/res-1.);
  histo,l0eff(where(l0eff<=l0lim)),ns;
  histo,l0g(where(l0g<=l0lim)),ns;
  histo,l0alt(where(l0alt<=l0lim)),ns;
  xytitles,"Outer scale (m)","Counts";
  lim = limits();
  plt,"A: Effective outer scale",lim(2)/2,lim(4)/2,tosys=1;
  plt,"B: Ground value ",lim(2)/2,lim(4)/2-lim(4)/8,tosys=1;
  plt,"C: Altitude value",lim(2)/2,lim(4)/2-lim(4)/4,tosys=1;
  pdf, "l0h_histo.pdf";
  w0 = where(l0eff<=l0lim);
  write,format="Average/1-sigma effective outer scale = %.3g m /%.3g m\n",avg(l0eff(w0)),(l0eff(w0))(rms);
  w0 = where(l0g<=l0lim);
  write,format="Average/1-sigma ground outer scale    = %.3g m /%.3g m\n",avg(l0g(w0)),(l0g(w0))(rms);
  w0 = where(l0alt<=l0lim);
  write,format="Average/1-sigma altitude outer scale  = %.3g m /%.3g m\n",avg(l0alt(w0)),(l0alt(w0))(rms);

  //Plotting altitude versus vh
  vlim=14;
  winkill,2;window,2,dpi=90,style="aanda.gs";
  plmk,alt(where(vh<vlim)),vh(where(vh<vlim)),marker=4,msize=0.1;
  range,-1,25;
  limits,-1,vlim;
  xytitles,"v(h) (m/s)","Altitude (km)";
  pdf, "vh_profile.pdf";

  //Plotting wind speed histograms
  winkill,6;window,6,dpi=90,style="aanda.gs";
  res = .4;
  ns = int(minmax(v0(where(v0<vlim)))(dif)/res-1.);
  histo,v0(where(v0<vlim)),ns;
  histo,v0g(where(v0g<vlim)),ns;
  histo,v0alt(where(v0alt<vlim)),ns;
  xytitles,"Wind speed (m/s)","Counts";
  lim = limits();
  plt,"A: Temporal coherence wind velocity",lim(2)/2,lim(4)/2,tosys=1;
  plt,"B: Ground value ",lim(2)/2,lim(4)/2-lim(4)/8,tosys=1;
  plt,"C: Altitude value",lim(2)/2,lim(4)/2-lim(4)/4,tosys=1;
  pdf, "vh_histo.pdf";
  
  w0 = where(v0g<vlim);
  write,format="Average/1-sigma global wind speed     = %.3g m/s /%.3g m/s\n",avg(v0(w0)),(v0(w0))(rms);
  w0 = where(v0g<vlim);
  write,format="Average/1-sigma ground wind speed     = %.3g m/s /%.3g m/s\n",avg(v0g(w0)),(v0g(w0))(rms);
  w0 = where(v0alt<vlim);
  write,format="Average/1-sigma altitude wind speed   = %.3g m/s /%.3g m/s\n",avg(v0alt(w0)),(v0alt(w0))(rms);
}

/*
 _____ ____  ____   ___  ____  
| ____|  _ \|  _ \ / _ \|  _ \ 
|  _| | |_) | |_) | | | | |_) |
| |___|  _ <|  _ <| |_| |  _ < 
|_____|_| \_\_| \_\\___/|_| \_\
                               
 ____  ____  _____    _    _  ______   _____        ___   _ 
| __ )|  _ \| ____|  / \  | |/ /  _ \ / _ \ \      / / \ | |
|  _ \| |_) |  _|   / _ \ | ' /| | | | | | \ \ /\ / /|  \| |
| |_) |  _ <| |___ / ___ \| . \| |_| | |_| |\ V  V / | |\  |
|____/|_| \_\_____/_/   \_\_|\_\____/ \___/  \_/\_/  |_| \_|
                                                            
*/
func plotScriptErrorBreakdown(Dir,scriptsuffix,path_out,redo=)
/* DOCUMENT eb = plotScriptErrorBreakdown("2013_09_17_onsky/","script293","results/");
 */
{
  local res;
  if(!redo)
    res = readfits(path_out + "scriptProcess_" + scriptsuffix +".fits",err=1);
  
  include, "PRAC_config.i",1;
  dataDirRoot = dataDir + Dir;
  //takes all slopestl files
  pathstl = listVersion(dataDirRoot,"fits","slopestl");
  pathstl = pathstl(where(strfind(scriptsuffix,pathstl)(2,) != -1));
  ntl = numberof(pathstl);
    
  if(is_void(res)){
    s_ir = s_tomo = s_noise = s_bw = s_fit = s_static = s_ncpa = s_ol = [];
    aomode = srsky = srmar = srpar = srborn  = mode = [];
    r0 = L0 = v =  h = ved1 = ved2 = ved3 = ved4 = ved5 = [];
    cnh2 = alt2 = l0h2 = vh2 = []; PSF = array(0.,70,70,ntl);
    res = array(pointer,24);
    
    for(j=1;j<=ntl;j++){
      timetl  = extractDate(pathstl(j));
      h       = grow(h,str2time(strpart(timetl,12:)));
      tmp     = PRAC_main(timetl,Dir = Dir,budgetonly=1);
      mode  = grow(mode,strchar(*tmp(1)));
      
      PSF(,,j)  = *tmp(5);
      
      p       = *tmp(4);
      s_ir  = grow(s_ir,p(1));
      s_fit  = grow(s_fit,p(2));
      s_bw  = grow(s_bw,p(3));
      s_tomo  = grow(s_tomo,p(4));
      s_alias  = grow(s_alias,p(5));
      s_noise  = grow(s_noise,p(6));
      s_ol  = grow(s_ol,p(7));
      s_static = grow(s_static,p(8));
      s_ncpa = grow(s_ncpa,p(9));
      srsky = grow(srsky,p(10));
      srmar = grow(srmar,p(11));
      srpar = grow(srpar,p(12));
      srborn = grow(srborn,p(13));
    
      r0 = grow(r0,(*tmp(2))(1));
      L0 = grow(L0,(*tmp(2))(2));
      v = grow(v,(*tmp(2))(3));

      ved1 = grow(ved1,(*tmp(6))(1));
      ved2 = grow(ved2,(*tmp(6))(2));
      ved3 = grow(ved3,(*tmp(6))(3));
      ved4 = grow(ved4,(*tmp(6))(4));
      ved5 = grow(ved5,(*tmp(6))(5));

      cnh2 = grow(cnh2,(*tmp(3))(,1));
      alt2 = grow(alt2,(*tmp(3))(,2));
      l0h2 = grow(l0h2,(*tmp(3))(,3));
      vh2  = grow(vh2,(*tmp(3))(,4));
      
    }
    res(1) = &mode;res(2) = &PSF;
    res(3) = &s_ir; res(4) = &s_tomo; res(5) = &s_alias; res(6) = &s_noise;
    res(7) = &s_bw;res(8) = &s_fit; res(9) = &s_ol;res(10) = &s_static;
    res(11) = &s_ncpa; res(12) = &srsky; res(13) = &srmar; res(14) = &srpar;res(15) = &srborn;
    res(16) = &r0;res(17) = &L0; res(18) = &v;res(19) = &h;
    res(20) = &[ved1,ved2,ved3,ved4,ved5];
    res(21) = &cnh2;
    res(22) = &alt2;
    res(23) = &l0h2;
    res(24) = &vh2;

    writefits,path_out + "scriptProcess_" + scriptsuffix +".fits",res;
  }

  //Displaying error breakdown
  modes = strchar(*res(1));
  sir_1  = (*res(3))(1::3);sir_2  = (*res(3))(2::3);sir_3  = (*res(3))(3::3);
  stomo_1  = (*res(4))(1::3);stomo_2  = (*res(4))(2::3);stomo_3  = (*res(4))(3::3);
  salias_1  = (*res(5))(1::3);salias_2  = (*res(5))(2::3);salias_3  = (*res(5))(3::3);
  snoise_1  = (*res(6))(1::3);snoise_2  = (*res(6))(2::3);snoise_3  = (*res(6))(3::3);
  sbw_1  = (*res(7))(1::3);sbw_2  = (*res(7))(2::3);sbw_3  = (*res(7))(3::3);
  sfit_1  = (*res(8))(1::3);sfit_2  = (*res(8))(2::3);sfit_3  = (*res(8))(3::3);
  sol_1  = (*res(9))(1::3);sol_2  = (*res(9))(2::3);sol_3  = (*res(9))(3::3);
  sstatic_1  = (*res(10))(1::3);sstatic_2  = (*res(10))(2::3);sstatic_3  = (*res(10))(3::3);
  sncpa_1  = (*res(11))(1::3);sncpa_2  = (*res(11))(2::3);sncpa_3  = (*res(11))(3::3);
  
  y1 = [sqrt(avg(sir_1^2)),sqrt(avg(stomo_1^2)),sqrt(avg(salias_1^2)),sqrt(avg(snoise_1^2)),sqrt(avg(sbw_1^2)),sqrt(avg(sfit_1^2)),sqrt(avg(sol_1^2)),sqrt(avg(sstatic_1^2)),sqrt(avg(sncpa_1^2))];
  y2 = [sqrt(avg(sir_2^2)),sqrt(avg(stomo_2^2)),sqrt(avg(salias_2^2)),sqrt(avg(snoise_2^2)),sqrt(avg(sbw_2^2)),sqrt(avg(sfit_2^2)),sqrt(avg(sol_2^2)),sqrt(avg(sstatic_2^2)),sqrt(avg(sncpa_2^2))];
  y3 = [sqrt(avg(sir_3^2)),sqrt(avg(stomo_3^2)),sqrt(avg(salias_3^2)),sqrt(avg(snoise_3^2)),sqrt(avg(sbw_3^2)),sqrt(avg(sfit_3^2)),sqrt(avg(sol_3^2)),sqrt(avg(sstatic_3^2)),sqrt(avg(sncpa_3^2))];
  
  labs = ["!s_IR","!s_tomo","!s_alias","!s_n","!s_bw","!s_fit","!s_ol","!s_stat","!s_ncpa"];
  winkill,0;window,0,style="aanda.gs",dpi=90;clr;
  plotsBarDiagram,y1,y2=y2,y3=y3,labs,col1=[char(241)],col2=[char(242)],col3=[char(243)],title=0;
  plt,modes(1),7,.8*max(y1(1),y2(1),y3(1)),tosys=1;
  plt,modes(2),7,.7*max(y1(1),y2(1),y3(1)),tosys=1,color=[128,128,128];
  plt,modes(3),7,.6*max(y1(1),y2(1),y3(1)),tosys=1,color=[64,64,64];
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/" + scriptsuffix +"/" +   scriptsuffix + "_budget.pdf";

  //Displaying averaged PSF

  PSF1 = ((*res(2))(,,1::3))(,,avg);
  PSF2 = ((*res(2))(,,2::3))(,,avg);
  PSF3 = ((*res(2))(,,3::3))(,,avg);
  a=3. ;
  cutmax = (max(max(PSF1),max(PSF2),max(PSF3)))^(1./a);
  box = 70; uz = 0.0297949; 

  window,1;clr;palette,"gray.gp";
  pli,(abs(PSF1))^(1./a),-box*uz/2,-box*uz/2,box*uz/2,box*uz/2,cmin = 0,cmax = cutmax;
  xytitles,"Arcsec","Arcsec";
  pltitle,"On-sky PSF in " + modes(1) + " with SR = " + var2str(arrondi(100*max(PSF1),1)) + "%";
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/" + scriptsuffix + "/"+ scriptsuffix +"_" + modes(1)+ "_psf.pdf";

  window,2;clr;palette,"gray.gp";
  pli,(abs(PSF2))^(1./a),-box*uz/2,-box*uz/2,box*uz/2,box*uz/2,cmin = 0,cmax = cutmax;
  xytitles,"Arcsec","Arcsec";
  pltitle,"On-sky PSF in " + modes(2) + " with SR = " + var2str(arrondi(100*max(PSF2),1)) + "%";
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/" + scriptsuffix + "/" + scriptsuffix +"_" + modes(2)+ "_psf.pdf";

  window,3;clr;palette,"gray.gp";
  pli,(abs(PSF3))^(1./a),-box*uz/2,-box*uz/2,box*uz/2,box*uz/2,cmin = 0,cmax = cutmax;
  xytitles,"Arcsec","Arcsec";
  pltitle,"On-sky PSF in " + modes(3) + " with SR = " + var2str(arrondi(100*max(PSF3),1)) + " %";
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/" + scriptsuffix + "/"+ scriptsuffix +"_" + modes(3)+ "_psf.pdf";

  //Displaying Strehl ratios versus time
  w1 = where(modes == modes(1))
  w2 = where(modes == modes(2));
  w3 = where(modes == modes(3));
  SRsky1 = (*res(12))(w1); SRsky2 = (*res(12))(w2); SRsky3 = (*res(12))(w3);
  SRmar1 = (*res(13))(w1); SRmar2 = (*res(13))(w2); SRmar3 = (*res(13))(w3);
  SRpar1 = (*res(14))(w1); SRpar2 = (*res(14))(w2); SRpar3 = (*res(14))(w3);
  SRborn1 = (*res(15))(w1); SRborn2 = (*res(15))(w2); SRborn3 = (*res(15))(w3);
  h1 = (*res(19))(w1); h2 = (*res(19))(w2); h3 = (*res(19))(w3);

  if(SRsky1(avg) >= 50){
    SRanal1 = SRmar1(sort(h1));
  }else if (50 > SRsky1(avg)  && SRsky1(avg) >= 15){
    SRanal1 = SRpar1(sort(h1));
  }else if (15 > SRsky1(avg)){
    SRanal1 = SRborn1(sort(h1));
  }

  if(SRsky2(avg) >= 50){
    SRanal2 = SRmar2(sort(h2));
  }else if (50 > SRsky2(avg) && SRsky2(avg) >= 15){
    SRanal2 = SRpar2(sort(h2));
  }else if (15 > SRsky2(avg)){
    SRanal2 = SRborn2(sort(h2));
  }

  if(SRsky3(avg) >= 50){
    SRanal3 = SRmar3(sort(h3));
  }else if (50 > SRsky3(avg) && SRsky3(avg) >= 15){
    SRanal3 = SRpar3(sort(h3));
  }else if (15 > SRsky3(avg)){
    SRanal3 = SRborn3(sort(h3));
  }

  SRanal1 = SRpar1(sort(h1));SRanal2 = SRpar2(sort(h2));SRanal3 = SRpar3(sort(h3));
  
  winkill,4;window,4,style="aanda.gs",dpi=90;clr;
  plmk, SRsky1(sort(h1)),h1(sort(h1)),msize=.4,marker=3,width=10; plg,SRsky1(sort(h1)),h1(sort(h1)),marks=0;
  plmk, SRanal1,h1(sort(h1)),msize=.4,marker=3; plg,SRanal1,h1(sort(h1)),marks=0,type=2;

  plmk, SRsky2(sort(h2)),h2(sort(h2)),msize=.4,marker=4,width=10; plg,SRsky2(sort(h2)),h2(sort(h2)),marks=0;
  plmk, SRanal2,h2(sort(h2)),msize=.4,marker=4; plg,SRanal2,h2(sort(h2)),marks=0,type=2;

  plmk, SRsky3(sort(h3)),h3(sort(h3)),msize=.4,marker=6,width=10; plg,SRsky3(sort(h3)),h3(sort(h3)),marks=0;
  plmk, SRanal3,h3(sort(h3)),msize=.4,marker=6; plg,SRanal3,h3(sort(h3)),marks=0,type=2;
  
  range,0,50;
  xytitles,"Local time (hour)","On-sky SR (%)";
  mm = min(grow(h1,h2,h3))
  plt,"Triangles: " + modes(1),mm,48,tosys=1;
  plt,"Circles: " + modes(2),mm,45,tosys=1;
  plt,"Crosses: " + modes(3),mm,42,tosys=1;
  plt,"Plain line: on-sky SR",mm+.05,48,tosys=1;
  plt,"Dash line: analytic SR",mm+.05,45,tosys=1;
  plt,"Script " + strpart(scriptsuffix,7:),mm+.05,42,tosys=1;
  limits,mm-0.01,mm+0.11;
  pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/"+ scriptsuffix + "/" + scriptsuffix +"_srvshour.pdf";
  //Displaying VED

   cnh = [median((*res(21))(1::5)),median((*res(21))(2::5)),median((*res(21))(3::5)),median((*res(21))(4::5)),median((*res(21))(5::5))];
   alt = [median((*res(22))(1::5)),median((*res(22))(2::5)),median((*res(22))(3::5)),median((*res(22))(4::5)),median((*res(22))(5::5))]/1000.;
    
  w1 = where(modes!="SCAO");
  if(is_array(w1)){
    m = modes(w1);
    w2 = w1(where(m==m(1)));
    w3 = w1(where(m==m(2)));

    ved1 = sqrt(((*res(20))(w2,)^2)(avg,));
    ved2 = sqrt(((*res(20))(w3,)^2)(avg,));
   
  
    winkill,5;window,5,style="aanda.gs",dpi=90;clr;
    thick=.3;
    for(j=1; j<=5; j++){
      d = 0.5*(0.5-thick);
      d = 2*thick;
      plfp, [char(242)],[0,1,1,0]*ved1(j),([1,1,1,1]*alt(j) -d/2) + [-1,-1,1,1]*thick,[4];
      plfp, [char(244)],[0,1,1,0]*ved2(j),([1,1,1,1]*alt(j) +d/2) + [-1,-1,1,1]*thick,[4];
    }

    xytitles,"Altitude (km)","VED (nm)";
    ml = max(grow(ved1,ved2));
    plt,m(1),1,ml*0.9,tosys=1,color=[128,128,128];
    plt,m(2),1,ml*0.8,tosys=1,color=[64,64,64];
    limits,min(alt)-1,max(max(alt)+1,15);
    range,0,ml*1.05;
    pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/"+ scriptsuffix + "/" + scriptsuffix + "_VED.pdf";
  
  
    //Displaying profiles
   
    
    winkill,6;window,6,style="aanda.gs",dpi=90;clr;
    wmoao = where(modes == "MOAO4L3N");
    suffcnh = readFitsKey(pathstl(wmoao(1)),"CN2H");
    ptrcnh  = restorefits("tomoparam",suffcnh);
    displayLayers,cnh,alt*1000,col=[char(242)],thick=.4,percent=1;
    displayLayers,-(*ptrcnh(2)),*ptrcnh(3),col=[char(241)],percent=1;
    limits,-100,100;
    range,-1,20;
    plt,"Retrieved profile\n on-sky",-90,15,tosys=1;
    plt,"Average profile\n during the script",20,15,tosys=1,color=[128,128,128];
    plt,"Script " + strpart(scriptsuffix,7:),-10,18,tosys=1;
    pdf,"/home/omartin/Documents/Articles/Canary_Phase_B/Figures/"+ scriptsuffix + "/" + scriptsuffix + "_profiles.pdf";
  }
  //gives statistics on turbulence
  //r0
  r0 = avg((*res(16))^(-5/3.))^(-3/5.);
  cnh2 = cnh*r0^(-5/3.)/sum(cnh);
  r0g = sum(cnh2(where(alt <= 0.1)))^(-3/5.);
  r0a = sum(cnh2(where(alt >0.1)))^(-3/5.);
  //seeing
  seeing = 0.103/r0;
  sg    = 0.103/r0g;
  sa    = 0.103/r0a;
  //Wind speed
  vh = [median((*res(24))(1::5)),median((*res(24))(2::5)),median((*res(24))(3::5)),median((*res(24))(4::5)),median((*res(24))(5::5))];
  v  = avg(vh); 
  vg = avg(vh(where(alt<=0.1)));
  va = avg(vh(where(alt>0.1)));

  //Outerscale
  alt = (*res(22));
  l0h = (*res(23));
  w = where(l0h <=100.);
  l0h = l0h(w);
  alt = alt(w);
  L0 =  avg(l0h^(5/3.))^(3/5.);
  L0g = avg((l0h(where(alt<=0.1)))^(5/3.))^(3/5.);
  L0a = avg((l0h(where(alt>0.1)))^(5/3.))^(3/5.);

  write,format="r0 (500nm)       = %.3g cm\n",100*r0;
  write,format="Ground r0        = %.3g cm\n",100*r0g;
  write,format="Altitude r0      = %.3g cm\n",100*r0a;
  write,format="Seeing           = %.3g ''\n",seeing;
  write,format="Ground seeing    = %.3g ''\n",sg;
  write,format="Altitude seeing  = %.3g ''\n",sa;
  write,format="Wind speed       = %.3g m/s\n",v;
  write,format="Ground windspeed = %.3g m/s\n",vg;
  write,format="Alt. windspeed   = %.3g m/s\n",va;
  write,format="L0               = %.3g m\n",L0;
  write,format="Ground L0        = %.3g m\n",L0g;
  write,format="Altitude L0      = %.3g m\n",L0a;
  
}




  
func pracAll(void)
/* DOCUMENT PRAC_all

 */
{
  
  

  tmpDir          = "/home/olivier/CANARY/Data/PHASE_B/tmp/stycomartin/";
  logFile         = tmpDir + "log.txt"; 
  dialogFile      = tmpDir + "dialogFile.txt"; 
  dataDir 	= "/home/olivier/CANARY/Data/PHASE_B/";

  Dir  = listFile("/home/olivier/CANARY/Data/PHASE_B/");
  Dir  = Dir(sort(Dir));
  Dir  = Dir(17:21);
  Dir  = "2012_10_04_onsky";
  ndir = numberof(Dir);
  
  goodDir = "results/pracResults/";
  if(!direxist(goodDir)) system,"mkdir " + goodDir;
  
  for(i=1;i<=ndir;i++){
    procDir = Dir(i);
    dataDirRoot = dataDir + procDir;
    //takes all slopestl files
    pathstl = listVersion(dataDirRoot,"fits","slopestl");
    if(is_array(pathstl)){
        //keeps only the script files
        w = where(strfind("script",pathstl)(0,) != -1);
        pathstl = pathstl(w);
        ntl = numberof(pathstl);
        //PRAC_res = array(pointer,32);
        //nl = 5;
        for(j=1;j<=ntl;j++){
          //extracts the date
          timetl = extractTime(pathstl(j));
          p = pracMain(timetl,Dir=procDir,verb=0,disp=0);
          /*if(is_pointer(p)){
            //Data identity
            for(k1=1;k1<=4;k1++){PRAC_res(k1) = & grow(*PRAC_res(k1),(strchar( (*p(1))(k1) )));}
            //global parameters
            for(k2=5;k2<=8;k2++){PRAC_res(k2) = & grow(*PRAC_res(k2),(*p(2))(k2-4));}
            //turbulence profile
            for(k3=9;k3<=11;k3++){PRAC_res(k3) = & grow(*PRAC_res(k3),(*p(3))((k3-9)*nl+1:(k3-8)*nl));}
            //tracking
            for(k8=12;k8<=14;k8++){PRAC_res(k8) = & grow(*PRAC_res(k8),(*p(4))(k8-11));}
            //error budget
            for(k4=15;k4<=21;k4++){PRAC_res(k4) = & grow(*PRAC_res(k4),(*p(5))(k4-14));}
            //Strehl ratios
            for(k5=22;k5<=31;k5++){PRAC_res(k5) = & grow(*PRAC_res(k5),(*p(6))(k5-21));}
            //Ensquared Energy
            PRAC_res(32) = &grow(*PRAC_res(32),*p(7));
            }*/
          write,format="%s","Job done: " + var2str(100.*j/ntl) + "%\r";
        }
        //writefits, goodDir + "PRACres_" + Dir(i) +".fits",PRAC_res;
    }
  }
}

/*
 ____  ____  _____     ____  
|  _ \/ ___||  ___|   |  _ \ 
| |_) \___ \| |_ _____| |_) |
|  __/ ___) |  _|_____|  _ < 
|_|   |____/|_|       |_| \_\
                             
*/


func plotSR(void)
{
  Dir = ["2013_09_17","2013_09_18","2013_09_16","2013_09_13","2013_09_15"];

  ndir = numberof(Dir);
  tmpDir          = "/mnt/CanaryData/tmp/stycomartin/";
  logFile         = tmpDir + "log.txt"; 
  dialogFile      = tmpDir + "dialogFile.txt"; 
  dataDir 	  = "/mnt/CanaryData/";
  

  SRf = SRn = SRsky = snr = aomode =  een=eef = [];
  
  for(i=1;i<=ndir;i++){
    //loading the results
    pf = readfits("results/PRACres_" + Dir(i)+".fits",err=1,verb=0);
    pn = readfits("results/PRACres_nofitl0_" + Dir(i)+".fits",err=1,verb=0);
    //AO modes
    aomode = grow(aomode,strchar(*pf(3)));
    //Strehl ratios
    SRsky = grow(SRsky,*pf(22));
    snr = grow(snr,*pf(23));
    SRf = grow(SRf,*pf(24));
    SRn = grow(SRn,*pn(24));
    eef = grow(eef,*pf(32));
    een = grow(een,*pn(32));
  }

  w = where(aomode == "MOAO" & SRsky>10);
  SRsky = SRsky(w);
  SRf = SRf(w);
  SRn = SRn(w);
  snr = snr(w);
  //Statistics
  write,format="Number of processed files: %d\n",numberof(SRsky);
  //Display
  winkill,0;
  window,0,style="aanda.gs",dpi=90;
  plmk,SRf,SRsky,msize=.3,width=4,marker=4,color="red";
  plg,[0,50],[0,50],marks=0,color="red";
  xytitles,"Sky strehl ratio","SR from PSR-R with L_0_(h) fitting";
  //Display
  winkill,1;
  window,1,style="aanda.gs",dpi=90;
  plmk,SRn,SRsky,msize=.3,width=4,marker=4,color="blue";
  plg,[0,50],[0,50],marks=0,color="black";
  xytitles,"Sky strehl ratio","SR from PSR-R with infinite altitude L_0_";
  //
  winkill,2;
  n = 200;
  window,2,style="aanda.gs",dpi=90;
  histo,eef,n,marks=0,color="red";
  histo,een,n,marks=0,color="blue";
  xytitles,"Accumulated error on the EE (%)","Counts";
  plt,"L_0_(h) fitting",10,60,color="red",tosys=1;
  plt,"no L_0_(h) fitting",10,50,color="blue",tosys=1;
  limits,-1,30;
  //
  winkill,3;
  window,3,style="aanda.gs",dpi=90;
  n = 15;
  af = 100*(SRf-SRsky)/SRsky;
  an = 100*(SRn-SRsky)/SRsky;
  wf= where(abs(af)<=100);
  wn= where(abs(an)<=100);
  histo,af(wf),n,color="red",marks=0;
  histo,an(wn),n,color="blue",marks=0;
  xytitles,"Error on the reconstructed SR (%)","Counts";
  plt,"L_0_(h) fitting",-100,60,color="red",tosys=1;
  plt,"no L_0_(h) fitting",-100,50,color="blue",tosys=1;
  avg(abs(af(wf)));
  avg(abs(an(wn)));

}


/*
 ____             __ _ _           
|  _ \ _ __ ___  / _(_) | ___  ___ 
| |_) | '__/ _ \| |_| | |/ _ \/ __|
|  __/| | | (_) |  _| | |  __/\__ \
|_|   |_|  \___/|_| |_|_|\___||___/
*/

func plotl0h(void,fitl0=)
/* plotl0h,fitl0=1;
   

 */
{
  Dir = ["2013_09_17","2013_09_18","2013_09_16","2013_09_13","2013_09_15"];

  ndir = numberof(Dir);
  tmpDir          = "/mnt/CanaryData/tmp/stycomartin/";
  logFile         = tmpDir + "log.txt"; 
  dialogFile      = tmpDir + "dialogFile.txt"; 
  dataDir 	  = "/mnt/CanaryData/";
  
  l0s = "nofitl0";
  if(fitl0) l0s = "" ;
  cnh = alt = l0h = cnhalt = [];
  
  for(i=1;i<=ndir;i++){
    procDir       = Dir(i) + "_onsky/";
    dataDirRoot   = dataDir + procDir;
    goodDir        = "profiles/" + procDir;
    //takes all slopestl files
    pathstl = listVersion(dataDirRoot,"fits","slopestl");
    //keeps only the script files
    w = where(strfind("script",pathstl)(0,) != -1);
    pathstl = pathstl(w);
    ntl = numberof(pathstl);

    for(j=1;j<=ntl;j++){
      //loading the prodiles
      h = extractDate(pathstl(j));
      p = readfits(goodDir + "ptrcn2h_"+h+"_nl_10"+l0s+".fits",err=1,verb=0);

      if(is_pointer(p)){
        //computing airmass
        object = strcase(1,readFitsKey(pathstl(j),"OBJECT"));
        tmp = givesCoordinatesFromObject(object);
        array_ra = str2flt(decoupe(tmp(1),' '));
        array_dec = str2flt(decoupe(tmp(2),' '));
        airm = airmassFromDate(extractDate(pathstl(j)),array_ra,array_dec);
        //unbiases from airmass
        tmpc = (*p(1))/airm;
        tmpa = (*p(2))/airm;
        //merges results
        cnh    = grow(cnh,tmpc);
        cnhalt = grow(cnhalt,sum(tmpc(2:)));
        alt    = grow(alt,tmpa);
        l0h    = grow(l0h,*p(3));
      }
    }
  }
  //Statistics
  w = where(alt > 0 & l0h <1e4);
  write,format="Number of processed files: %d\n",numberof(l0h)/10;
  write,format="Altitude layers with L_0<1e4m: %.3g%s\n",100.*numberof(w)/numberof(l0h(where(alt>0))),"%";
  //Display
  winkill,0;
  window,0,style="aanda.gs",dpi=90;
  plmk,alt/1000.,log10(l0h),msize=.1;
  xytitles,"log_10_(L_0_ (m))","Altitude (km)";
  range,-1,20;

  La = l0h(w);
  cnha = cnh(w);
  winkill,1;
  window,1,style="aanda.gs",dpi=90;
  plmk,log10(La),0.103/cnha^(-3/5.), msize=.1;
  xytitles,"Altitude seeing (arcsec)","log_10_(Altitude L_0_<1e4m (m))";
  limits,-.1,0.5;
  range,-2,5;
  w = where(alt > 0 & l0h <=.5e4 & 0.103/cnh^(-3/5.) > .05);
  write,format="Altitude layers with L_0<1e4m and a significant strength: %.3g%s\n",100.*numberof(w)/numberof(l0h(where(alt>0))),"%";

  winkill,2;
  window,2,style="aanda.gs",dpi=90;
  histo, 0.103/cnh(where(alt==0))^(-3/5.),50,marks=0,color="red";
  histo, 0.103/cnhalt^(-3/5.),100,marks=0,color="blue";
  limits,0,2;
  xytitles,"Seeing (arcsec)","Counts";

  winkill,3;
  window,3,style="aanda.gs",dpi=90;
  histo, l0h(where(alt == 0 & l0h<50)),100,marks=0,color="red";
  xytitles,"Ground Outer scale (m)","Counts";
}

func givesCoordinatesFromObject(object)
{
  object = strcase(1,object);

  if(object=="A47" | object=="AST47"){
    alpha  = "21 12 3.0";
    delta = "38 36 45.0";
  }
  if(object=="A453" | object=="AST453"){
    alpha  = "+1 36 44.400";
    delta = "47 21 23.0";
  }
  if(object=="A396" | object=="AST396"){
    alpha  = "6 48 15.0";
    delta = "41 05 44.0";
  }
  if(object=="A53" | object=="AST53"){
    alpha  = "+23 24 30.0";
    delta = "40 53 55.0";
  }
 if(object=="A110" | object=="AST110"){
    alpha  = "+19 31 11";
    delta = "34 56 02";
  }
 if(object=="AS2" | object=="A32" | object=="AST32"){
    alpha  = "+18 27 7.0";
    delta = "26 52 34.0";
  }
 if(object=="A34" | object=="AST34"){
    alpha  = "+18 51 30.0";
    delta = "10 20 3.0";
  }
 if(object=="A51" | object=="AST51"){
    alpha  = "+22 48 5";
    delta = "39 17 3.0";
  }

 if(object=="A10" | object=="AST10"){
    alpha  = "+5 52 17.0";
    delta = "32 34 36.0";
  }
 
 if(object=="AT1" | object=="ASTT1"){
    alpha  = "+13 41 38.2";
    delta = "+7 36 21.0";
  }

 if(object=="A12" | object=="AST12"){
    alpha  = "+6 1 9.0";
    delta = "+23 20 29.0";
  }
 
 
  return [alpha,delta];

}

func airmassFromDate(date,alphaJ2000,deltaJ2000)
/* DOCUMENT
 */
{
  //extract date from path
  Y = str2int(strpart(date,1:4));
  M = str2int(strpart(date,6:7));
  D = str2int(strpart(date,9:10));
  H = str2int(strpart(date,12:13))-1.;//correction temps local vers temps UTC
  Mm = str2int(strpart(date,15:16));
  S = str2int(strpart(date,18:19));
  //defines coordinates
  alpha = alpha2deg(alphaJ2000(1),alphaJ2000(2),alphaJ2000(3));
  delta = dec2deg(deltaJ2000(1),deltaJ2000(2),deltaJ2000(3));
  
  return airmassFromCoordinates(Y, M, D, H, Mm, S, alpha, delta);
}
func alpha2deg(h, m, s) {
        return dec2deg(h, m, s)*15;
}

func dec2deg(h, m, s) {
        return h+m/60.+s/3600.;
}
func airmassFromCoordinates(Y, M, D, H, Mm, S, alpha, delta, latitude, longitude)
/* DOCUMENT airmass(Y, M, D, H, Mm, S, alpha, delta, latitude, longitude)
     
   SEE ALSO:
 */
{
        LATITUDE = 28.7567;
        LONGITUDE = -17.8917;

        if (is_void(latitude)) latitude=LATITUDE;
        if (is_void(longitude)) longitude=LONGITUDE;

        alt = AAAltitude(Y, M, D, H, Mm, S, alpha, delta, latitude, longitude)(1);
        return (1/cos(pi/2.-alt*pi/180));
}
func AAAltitude(Y, M, D, H, Mm, S, alpha, delta, latitude, longitude) {
        LATITUDE = 28.7567;
        LONGITUDE = -17.8917;

        if (is_void(latitude)) latitude=LATITUDE;
        if (is_void(longitude)) longitude=LONGITUDE;


        lst = LST(Y, M, D, H, Mm, S, longitude) * 15.;
        H = lst - alpha;
        r = pi/180;
        sina = sin(latitude*r)*sin(delta*r)+cos(latitude*r)*cos(delta*r)*cos(H*r);
        alt = asin(sina) / r;
        return alt;
}
func LST(Y, M, D, H, Mm, S, longitude) {

        LONGITUDE = -17.8917;
        if (is_void(longitude)) longitude=LONGITUDE;

        JJ = JD(Y, M, D, H, Mm, S);
        D = JJ - 2451545.;
        GMST = 18.697374558 + 24.06570982441908 * D;
        lst = GMST+longitude/15.;
        lst%=24;
        return lst;
}

func JD(Y, M, D, H, Mm, S) {
        h=(H-12)/24. + Mm/24./60 + S/24./3600;
        J=(1461 * (Y + 4800 + (M - 14)/12))/4 +(367 * (M - 2 - 12 * ((M - 14)/12)))/12 - (3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075;
        return double(J) + h;
}

/*
 ____  _        _   _                        _ _         
/ ___|| |_ __ _| |_(_) ___  _ __   __ _ _ __(_) |_ _   _ 
\___ \| __/ _` | __| |/ _ \| '_ \ / _` | '__| | __| | | |
 ___) | || (_| | |_| | (_) | | | | (_| | |  | | |_| |_| |
|____/ \__\__,_|\__|_|\___/|_| |_|\__,_|_|  |_|\__|\__, |
                                                   |___/ 
*/


func plotcnhvstime(Dir)
/* plotcnhvstime,"2013_09_17"
   

 */
{

  tmpDir          = "/mnt/CanaryData/tmp/stycomartin/";
  logFile         = tmpDir + "log.txt"; 
  dialogFile      = tmpDir + "dialogFile.txt"; 
  dataDir 	  = "/mnt/CanaryData/";
  procDir         = Dir + "_onsky/";
  dataDirRoot     = dataDir + procDir;
  goodDir         = "profiles/" + procDir;

  cnh0 = cnha = L0 = t = [];
  
  //takes all slopestl files
  pathstl = listVersion(dataDirRoot,"fits","slopestl");
  //keeps only the script files
  w = where(strfind("script",pathstl)(0,) != -1);
  pathstl = pathstl(w);
  ntl = numberof(pathstl);
  
  for(j=1;j<=ntl;j++){
    //loading the prodiles
    h = extractDate(pathstl(j));
    p = readfits(goodDir + "ptrcn2h_"+h+"_nl_10nofitl0.fits",err=1,verb=0);

    if(is_pointer(p)){
      //computing airmass
      object = strcase(1,readFitsKey(pathstl(j),"OBJECT"));
      tmp = givesCoordinatesFromObject(object);
      array_ra = str2flt(decoupe(tmp(1),' '));
      array_dec = str2flt(decoupe(tmp(2),' '));
      airm = airmassFromDate(extractDate(pathstl(j)),array_ra,array_dec);
      //unbiases from airmass
      tmpc = (*p(1))/airm;
      //tmpa = (*p(2))/airm;
      //merges results
      cnh0 = grow(cnh0,(tmpc)(1));
      cnha = grow(cnha,sum((tmpc)(2:)));
      L0 = grow(L0,(*p(3))(1));
      t   = grow(t,str2time(strpart(h,12:)));
    }
  }
  
  //Display
  winkill,0;
  window,0,style="aanda.gs",dpi=90;
  plmk,0.103/cnh0^(-3/5.),t,msize=.1;
  xytitles,"Local time (hour)","Ground Seeing (arcsec)";
  range,0,1;
  
  winkill,1;
  window,1,style="aanda.gs",dpi=90;
  plmk,0.103/cnha^(-3/5.),t,msize=.1;
  xytitles,"Local time (hour)","Altitude Seeing (arcsec)";
  range,0,1;
  
  winkill,2;
  window,2,style="aanda.gs",dpi=90;
  plmk,L0,t,msize=.1;
  xytitles,"Local time (hour)","Ground L_0_";
  range,0,30;
}

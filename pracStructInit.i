
struct magnitude_struct{
  double H;
  double K;
  double Kp;
  double J;
  double V;
  double R;
  double B;

};
struct asterism_struct{
  string name;
  string rightAscension;
  string declinaison;
  double nOffAxisSources;
  double photometry;
  magnitude_struct magnitude(nngs);
};

struct atm_struct{
  double  lambda;
  double  r0;
  double  L0;
  double  v;
  double  dir;
  int     nLayers;
  double  cnh(nl);
  double  altitude(nl);
  double  l0h(nl);
  double  vh(nl);
  double  dirh(nl);
};

struct telescope_struct{
  double diam;
  double obs;
  double airmass;
  double zenith;
  double azimuth;
  double plateScale;
  long   nPixels;
  double fourierPixSize;
  double pixSize;
  double foV;
  double pitch;
  long   pitchSize;
  double lambdaPerFov;
  double dPerFov;
  double airyPeak;
};

struct rtc_struct{
  int      nFrames; 
  double   loopGain;
  double   delay;
  double   frameDelay;
  double   Fe;
  double   BP;
  string   obsMode;
  string   recType;
  string   aoMode;

  int      nWfs;
  int      nLgsWfs;
  int      nNgsWfs;
  int      nSlopes;
  int      its; // index of truth sensor
  pointer  ptrListOffAxisNgs;     
  pointer  ptrListLgs;

  pointer  R;
  pointer  skyProfile;
  pointer  mi;
  pointer  mc;
  
  pointer  slopes_res;
  pointer  slopes_dis;
  pointer  volts;
};

struct wfs_struct{
  double x; // This WFS guide star x-position in arcsec. Optional [0,0]
  double y; // This WFS guide star y-position in arcsec. Optional [0,0]
  double lgsH; // This WFS GS altitude in meter. 0 for infinity. Specified at zenith. Optional [0]
  double lgsdH;// This WFS GS depth in meter (e.g. Na layer thickness).
  double lgsAngle;
  double lgsRadius;
  long   nLenslet;
  int    nValidSubAp;
  int    nValidMeas;
  double pixSize;        // Subaperture pixel size in arcsec. Required [none]
  int    nPixels;        // Final # of pixels in subaperture. Required [none]
  int    type; // 0 if unused. 1 if NGS (HO + TT) , 2 if LGS (HO), 3 (TT only)
  long   sym;
  double unit;
  int    flip;
};

struct dm_struct {
  long    nActu;
  long    nValidActu;
  double  pitch;
  pointer csX;
  pointer csY;
  pointer csI;  // indices of actuators = [3,4,5,6, 10,11,12,13,14,15, 17 ... ] in a Ndiam X Ndiam image
  double  x0;
  double  couplage;
};

struct cam_struct{
  string  name;
  pointer psfRaw;
  pointer bg;
  pointer psfCal;
  long    nxPixels;   // dimension along X
  long    nyPixels;   // dimension along Y
  long    xPsf, yPsf;
  double  lambda;
  double  pixSize;
  double  dPixSizeOnLambda;
  pointer deadPixIndex;
  int     nPixelsCropped;
  double  exposureTime;
};

struct sys_struct{
  double  tracking(3);
  double  xshift(nwfs);
  double  yshift(nwfs);
  double  magnification(nwfs);
  double  theta(nwfs);
  double  centroidGain(nwfs);
  pointer slopesToZernikeMatrix;
  pointer zernikeToSlopesMatrix;
};

struct learn_struct{
  int     ttr;
  int     diagonal;
  pointer transformationMatrix;
  int     runLearnTwoSteps;
  int     nl;
  double  cnh(nl);
  double  altitude(nl);
  double  l0h(nl);
  double  tracking(3);
  double  xshift(nwfs);
  double  yshift(nwfs);
  double  magnification(nwfs);
  double  theta(nwfs);
  double  centroidGain(nwfs)
};

struct learnfit_struct{
  int     cnh(nl);
  int     altitude(nl);
  int     l0h(nl);
  int     tracking(3);
  int     xshift(nwfs);
  int     yshift(nwfs);
  int     magnification(nwfs);
  int     theta(nwfs);
  int     centroidGain(nwfs);
};

struct covmatrix_struct{
  pointer slopes;
  pointer learn;
  pointer aliasing;
  pointer parallel;
  pointer noise;
  pointer tracking;
};

struct budget_struct{
  double  res;
  double  fit;
  double  bw;
  double  tomo;
  double  alias;
  double  noise;
  double  static;
  double  ved(nl);
  double  ol;
  double  ncpa;
  double  vib
  double  SRsky;
  double  SRres;
  double  SRmar;
  double  SRpar;
  double  SRborn;
  double  SRncpa;
};

struct otf_struct{
  pointer tel;
  pointer cpa;
  pointer fit;
  pointer bw;
  pointer tomo;
  pointer ncpa;
  pointer static;
  pointer ts;
  pointer ts_cropped;
  pointer sky;
  pointer res;
};

struct psf_struct{
  pointer sky;
  pointer res;
  pointer diff;
  pointer ncpa;
  pointer EE_res;
  pointer EE_sky;
  double  SR_res;
  double  SR_sky;
  double  SR_tomo;
  double  SR_bw;
  double  SR_fit;
  double  SR_stats;
  double  SR_ncpa;
  double  diffSum;
  double  diffRms;
  double  FWHM_res;
  double  FWHM_sky;
  
}

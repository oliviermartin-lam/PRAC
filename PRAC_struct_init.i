
struct telescope_struct{
  double diam;
  double obs;
  double airm;
  double zen;
  string object;
};

struct wfs_struct{
  int      nssp;
  long     sX;           // linear number of subapertures (7)
  long     sY;
  int      type; // 0 if unused. 1 if NGS (HO + TT) , 2 if LGS (HO), 3 (TT only)
  double   pixarc
  double   lgsH;
  double   lgsdH;   // thickness of sodium layer
  double   lgsAngle;
  double   lgsRadius;
  double   x; // x pos in ''
  double   y; // y pos in ''
  long     sym;
  double   unit;
  pointer  mrz;
  pointer  zmi;
  pointer  filterouttiltMatrix;  // matrix that filters tilt out from slopes measurements
  pointer  keeponlytiltMatrix;  // matrix that just keep tiptilt from slopes measurements
  int      flip;
};

struct dm_struct {
  long    nactu;
  long    ndiam;
  pointer csX;
  pointer csY;
  pointer csI;  // indices of actuators = [3,4,5,6, 10,11,12,13,14,15, 17 ... ] in a Ndiam X Ndiam image
  double x0;
};

struct rtc_struct{
  double   gain;
  double   delay;
  double   Fe;
  double   BP;
  string   obsmode;
  string   rectype;
  string   aomode;
  pointer  R;
  pointer  ptrlistOffAxisNgs;     
  pointer  ptrlistLgs;
  pointer  MI;
  pointer  MC;
  pointer  volts;
};

struct ircam_struct {
  string      camName;
  pointer     psfraw;
  pointer     bg;
  pointer     psfcal;
  long        dx;   // dimension along X
  long        dy;   // dimension along Y
  long        xPsf, yPsf;
  double      lambda_ir;
  double      pixSize;
  double      uz;
  double      uld;
  pointer     deadPixIndex;
  int         npix_sr;
  
};

struct fourier_struct {
  long    npix;
  double  uk;
  double  ud;
  double  champ_Dphi;
  double  dactu;
  long    dactupix;
  double  uz;
  double  uld;
  pointer k;
  pointer kx;
  pointer ky;
};

struct turbu_struct{
  double  r0ir;
  double  r0vis;
  double  L0;
  double  v;
  double  dir;
  pointer varNoise;
};
struct learn_struct{
  int     ttr;
  int     diagonal;
  pointer transformationMatrix;
  int     runLearnTwoSteps;
  int     nl;
  double  cnh(NL_DEF);
  double  altitude(NL_DEF);
  double  l0h(NL_DEF);
  double  vh(NL_DEF);
  double  dirh(NL_DEF);
  double  tracking(3);
  double  xshift(NBWFS_DEF);
  double  yshift(NBWFS_DEF);
};
struct learnfit_struct{
  int     cnh(NL_DEF);
  int     altitude(NL_DEF);
  int     l0h(NL_DEF);
  int     tracking(3);
  int     xshift(NBWFS_DEF);
  int     yshift(NBWFS_DEF);
};
struct covmatrix_struct{
  pointer covMeas;
  pointer covLearn;
  pointer covAlias;
  pointer covPara;
  pointer covNoise;
  pointer covTracking;
};
struct budget_struct{
  double  res;
  double  fit;
  double  bw;
  double  tomo;
  double  alias;
  double  noise;
  double  static;
  double  ved(NL_DEF);
  double  ol;
  double  ncpa;
  double  vib
  double  SRsky;
  double  SRres;
  double  SRmar;
  double  SRpar;
  double  SRborn;
};
struct uncertainties_struct{
  double     r0;
  double     L0;
  double     v;
  double     dir;
  double     tracking(3);
  double     cnh(NL_DEF);
  double     l0h(NL_DEF);
  double     altitude(NL_DEF);
  double     vh(NL_DEF);
  double     dirh(NL_DEF);
  double     srsky;
};
struct data_struct{
  int     nwfs;
  double  lambda_vis;
  int     its;                // index of truth sensor
  int     nframes;
  pointer slopes_res;
  pointer slopes_dis;
  int     Nslopes;
  double  plscale;
  
  wfs_struct wfs(NBWFS_DEF);
  telescope_struct tel;
  fourier_struct fourier;
  dm_struct dm;
  ircam_struct camir;
  rtc_struct rtc;
  turbu_struct turbu;
  learn_struct learn;
  learnfit_struct learn_fit;
  budget_struct budget;
  uncertainties_struct uncertainties;
  covmatrix_struct covmatrix;
};

data = data_struct();

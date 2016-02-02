/*
 ___ _   _ ____ _____  _    _   _ _____  _    _   _ _____ ___  _   _ ____  
|_ _| \ | / ___|_   _|/ \  | \ | |_   _|/ \  | \ | | ____/ _ \| | | / ___| 
 | ||  \| \___ \ | | / _ \ |  \| | | | / _ \ |  \| |  _|| | | | | | \___ \ 
 | || |\  |___) || |/ ___ \| |\  | | |/ ___ \| |\  | |__| |_| | |_| |___) |
|___|_| \_|____/ |_/_/   \_\_| \_| |_/_/   \_\_| \_|_____\___/ \___/|____/ 
                                                                           
 ____  ____  _____     ____  
|  _ \/ ___||  ___|   |  _ \ 
| |_) \___ \| |_ _____| |_) |
|  __/ ___) |  _|_____|  _ < 
|_|   |____/|_|       |_| \_\

*/

func mmseParallelModes(slopes)
/* DOCUMENT

 */

{

  g = slopes - slopes(,avg);

  
  //empirical covariance matrix of the residue
  C_all = g(,+) * g(,+)/dimsof(g)(0);
  
  //inversion
  C_all_m = invgen(C_all,0,cond=30.);

  //closed-loop noise matrix
  varNoise  = getNoiseVar(g);//in arcsec^2
  Cnn       = 0*C_all;
  takesDiag(Cnn) = varNoise;
  
  //aliasing matrix
  if(rtc.obsMode == "MOAO"){

    Crr    = *covMatrix.aliasing;
    Crr    = computeCovSlopesError(Crr,*rtc.R);
    
  } else if(rtc.obsMode == "SCAO"){

    Crr = (*covMatrix.aliasing)(slrange(rtc.its),slrange(rtc.its));
  }

  C_ee = C_all - Cnn  + Crr;

  //MMSE reconstruction
  Rmmse  = C_ee(,+) * C_all_m(+,);
  gpara    = Rmmse(,+) * g(+,); // in arcsec

  return gpara;
}

/*
__     _______ ____      _    _   _   ____  ____  _____     ____  
\ \   / / ____|  _ \    / \  | \ | | |  _ \/ ___||  ___|   |  _ \ 
 \ \ / /|  _| | |_) |  / _ \ |  \| | | |_) \___ \| |_ _____| |_) |
  \ V / | |___|  _ <  / ___ \| |\  | |  __/ ___) |  _|_____|  _ < 
   \_/  |_____|_| \_\/_/   \_\_| \_| |_|   |____/|_|       |_| \_\

*/

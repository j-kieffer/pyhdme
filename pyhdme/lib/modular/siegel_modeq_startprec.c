
#include "modular.h"

slong siegel_modeq_startprec(fmpz* I, slong ell)
{
  slong weights[4] = IGUSA_WEIGHTS;
  slong h = cov_height(I, 4, weights);
  slong d = siegel_nb_cosets(ell);
  
  slong res = SIEGEL_START_PREC_MUL * (n_pow(ell, 3) * n_clog(ell, 2) + d*h)
    + SIEGEL_START_PREC_ADD;
  return 100 * (res/100 + 1);
}

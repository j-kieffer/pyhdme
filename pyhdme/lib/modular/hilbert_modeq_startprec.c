
#include "modular.h"

slong hilbert_modeq_startprec(fmpz* I, slong ell, slong delta)
{
  slong weights[4] = IGUSA_WEIGHTS;
  slong h = cov_height(I, 4, weights);
  slong d = hilbert_nb_cosets(ell, delta);
  
  slong res = HILBERT_START_PREC_MUL * ell * n_clog(ell, 2)
    + 20 * d * h
    + HILBERT_START_PREC_ADD;
  return 100 * (res/100 + 1);
}

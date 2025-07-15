
#include "hecke.h"

slong siegel_nb_T1_cosets(slong ell)
{
  return ell + n_pow(ell, 2) + n_pow(ell, 3) + n_pow(ell, 4);
}


#include "hecke.h"

slong siegel_nb_cosets(slong ell)
{
  return (n_pow(ell,4) - 1) / (ell - 1);
}

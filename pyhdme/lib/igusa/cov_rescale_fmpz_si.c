
#include "igusa.h"

void cov_rescale_fmpz_si(fmpz* I, fmpz* S, slong scal, slong nb, slong* weights)
{
  fmpz_t f;

  fmpz_init(f);
  fmpz_set_si(f, scal);
  cov_rescale_fmpz(I, S, f, nb, weights);
  
  fmpz_clear(f);
}


#include "igusa.h"

int cov_divisible_fmpz(fmpz* I, const fmpz_t scal, slong nb, slong* weights)
{
  fmpz_t f;
  int r = 1;
  slong j;

  fmpz_init(f);

  for (j = 0; j < nb; j++)
    {
      fmpz_pow_ui(f, scal, weights[j]);
      r = r && fmpz_divisible(&I[j], f);
    }

  fmpz_clear(f);
  return r;
}

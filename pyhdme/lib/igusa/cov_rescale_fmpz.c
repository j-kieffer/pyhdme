
#include "igusa.h"

void cov_rescale_fmpz(fmpz* I, fmpz* S, const fmpz_t scal, slong nb, slong* weights)
{
  fmpz* res;
  fmpz_t f;
  slong j;

  res = _fmpz_vec_init(nb);
  fmpz_init(f);

  for (j = 0; j < nb; j++)
    {
      fmpz_pow_ui(f, scal, weights[j]);
      fmpz_mul(&res[j], &S[j], f);
    }

  _fmpz_vec_set(I, res, nb);  
  _fmpz_vec_clear(res, nb);
  fmpz_clear(f);
}

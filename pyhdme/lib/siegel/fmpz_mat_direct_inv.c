
#include "siegel.h"

void fmpz_mat_direct_inv(fmpz_mat_t minv, const fmpz_mat_t m)
{
  fmpz_mat_t res;
  fmpz_t den;
  
  fmpz_init(den);
  fmpz_mat_init(res, fmpz_mat_nrows(m), fmpz_mat_ncols(m));

  fmpz_mat_inv(res, den, m);
  if (!fmpz_is_one(den))
    {
      fmpz_mat_neg(res, res);
      fmpz_neg(den, den);
    }
  if (!fmpz_is_one(den))
    {
      flint_printf("(fmpz_mat_direct_inv) Error: impossible inverse\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_mat_set(minv, res);
  
  fmpz_mat_clear(res);
  fmpz_clear(den);
}

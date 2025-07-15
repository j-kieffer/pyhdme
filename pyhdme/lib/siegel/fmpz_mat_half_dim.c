
#include "siegel.h"

slong fmpz_mat_half_dim(const fmpz_mat_t m)
{
  slong g = fmpz_mat_nrows(m)/2;
  if (fmpz_mat_nrows(m) != 2*g
      || fmpz_mat_ncols(m) != 2*g)
    {
      flint_printf("(fmpz_mat_half_dim) Wrong matrix dimensions\n");
      fflush(stdout);
      flint_abort();
    }
  return g;
}

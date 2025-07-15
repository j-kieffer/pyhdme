
#include "siegel.h"

void fmpz_mat_diagonal_symplectic(fmpz_mat_t m, const fmpz_mat_t u)
{
  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t d, zero;
  fmpz_t den;

  fmpz_mat_init(d, g, g);
  fmpz_mat_init(zero, g, g);
  fmpz_init(den);

  fmpz_mat_inv(d, den, u);
  fmpz_mat_transpose(d, d);
  if (!fmpz_is_one(den))
    {
      fmpz_neg(den, den);
      fmpz_mat_neg(d, d);
    }
  if (!fmpz_is_one(den))
    {
      flint_fprintf(stderr, "(fmpz_mat_diagonal_symplectic) Error: non-invertible matrix\n");
      flint_abort();
    }

  fmpz_mat_set_abcd(m, u, zero, zero, d);

  fmpz_mat_clear(d);
  fmpz_mat_clear(zero);
  fmpz_clear(den);
}

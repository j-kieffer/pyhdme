
#include "siegel.h"

int fmpz_mat_is_J(const fmpz_mat_t m)
{
  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t x;
  int res = 1;

  fmpz_mat_init(x, g, g);
  
  fmpz_mat_get_a(x, m);
  if (!fmpz_mat_is_zero(x)) res = 0;
  fmpz_mat_get_b(x, m);
  if (!fmpz_mat_is_one(x)) res = 0;
  fmpz_mat_get_c(x, m);
  fmpz_mat_neg(x, x);
  if (!fmpz_mat_is_one(x)) res = 0;
  fmpz_mat_get_d(x, m);
  if (!fmpz_mat_is_zero(x)) res = 0;

  fmpz_mat_clear(x);
  return res;
}

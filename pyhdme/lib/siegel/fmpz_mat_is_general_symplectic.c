
#include "siegel.h"

int fmpz_mat_is_general_symplectic(const fmpz_mat_t m)
{
  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t r, J;
  int res;

  fmpz_mat_init(r, 2*g, 2*g);
  fmpz_mat_init(J, 2*g, 2*g);
  fmpz_mat_J(J);
  
  fmpz_mat_transpose(r, m);
  fmpz_mat_mul(r, r, J);
  fmpz_mat_mul(r, r, m);
  fmpz_mat_mul(r, r, J);

  res = fmpz_mat_is_scalar(r)
    && !fmpz_is_zero(fmpz_mat_entry(r, 0, 0));

  fmpz_mat_clear(r);
  fmpz_mat_clear(J);
  return res;
}

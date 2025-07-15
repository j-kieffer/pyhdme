
#include "hilbert.h"

int hilbert_is_totally_positive(const fmpz_poly_t x, slong delta)
{
  fmpz_t lhs, rhs;
  int res = 1;
  
  fmpz_init(lhs);
  fmpz_init(rhs);

  fmpz_poly_get_coeff_fmpz(lhs, x, 0);
  fmpz_mul_si(lhs, lhs, 2);
  fmpz_poly_get_coeff_fmpz(rhs, x, 1);
  if (delta % 2 == 1)
    {
      fmpz_add(lhs, lhs, rhs);
    }
  if (fmpz_cmp_si(lhs, 0) <= 0)
    {
      res = 0;
    }
  
  fmpz_mul(lhs, lhs, lhs);
  fmpz_mul(rhs, rhs, rhs);
  fmpz_mul_si(rhs, rhs, delta);
  if (fmpz_cmp(lhs, rhs) <= 0)
    {
      res = 0;
    }

  fmpz_clear(lhs);
  fmpz_clear(rhs);
  return res;
}

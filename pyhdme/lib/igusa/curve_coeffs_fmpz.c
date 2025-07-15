
#include "igusa.h"

void curve_coeffs_fmpz(fmpz* ai, const fmpz_poly_t crv)
{
  slong k;
  fmpz_t coeff;

  fmpz_init(coeff);
  for (k = 0; k < 7; k++)
    {
      fmpz_poly_get_coeff_fmpz(coeff, crv, k);
      fmpz_set(&ai[k], coeff);
    }
  fmpz_clear(coeff);
}

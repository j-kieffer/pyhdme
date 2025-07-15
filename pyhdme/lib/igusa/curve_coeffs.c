
#include "igusa.h"

void curve_coeffs(acb_ptr ai, const acb_poly_t crv)
{
  slong k;
  acb_t coeff;

  acb_init(coeff);
  for (k = 0; k < 7; k++)
    {
      acb_poly_get_coeff_acb(coeff, crv, k);
      acb_set(&ai[k], coeff);
    }
  acb_clear(coeff);
}


#include "hilbert.h"

void hilbert_sigma1(acb_t z, const fmpz_poly_t x, slong delta, slong prec)
{
  acb_t res;
  fmpz_t coeff;
  
  acb_init(res);
  fmpz_init(coeff);

  fmpz_poly_get_coeff_fmpz(coeff, x, 0);
  acb_set_fmpz(z, coeff);

  acb_zero(res);
  arb_sqrt_ui(acb_realref(res), delta, prec);
  if (delta % 2 == 1)
    {
      acb_add_si(res, res, 1, prec);
    }
  acb_div_si(res, res, 2, prec);

  fmpz_poly_get_coeff_fmpz(coeff, x, 1);
  acb_mul_fmpz(res, res, coeff, prec);
	       
  acb_add(z, z, res, prec);
  
  acb_clear(res);
  fmpz_clear(coeff);
}

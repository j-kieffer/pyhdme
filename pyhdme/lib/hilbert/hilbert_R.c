
#include "hilbert.h"

void hilbert_R(acb_mat_t R, slong delta, slong prec)
{
  fmpz_poly_t x;
  fmpz_poly_init(x);

  fmpz_poly_one(x);
  hilbert_sigma1(acb_mat_entry(R, 0, 0), x, delta, prec);
  hilbert_sigma2(acb_mat_entry(R, 1, 0), x, delta, prec);

  fmpz_poly_zero(x);
  fmpz_poly_set_coeff_si(x, 1, 1);
  hilbert_sigma1(acb_mat_entry(R, 0, 1), x, delta, prec);
  hilbert_sigma2(acb_mat_entry(R, 1, 1), x, delta, prec);
  
  fmpz_poly_clear(x);
}

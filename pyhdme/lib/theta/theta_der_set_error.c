
#include "theta.h"

int theta_der_set_error(mag_t error, const acb_mat_t tau, slong prec)
{
  arb_mat_t im;
  arb_t lambda;
  arb_t one;
  arf_t ubound;
  fmpz_t a, b;
  int res;

  arb_mat_init(im, 2, 2);
  arb_init(lambda);
  arb_init(one);
  arf_init(ubound);
  fmpz_init(a);
  fmpz_init(b);

  acb_mat_get_imag(im, tau);
  arb_mat_lambda(lambda, im, prec);
  arb_inv(lambda, lambda, prec);
  arb_one(one);
  arb_max(lambda, lambda, one, prec);  
  arb_get_ubound_arf(ubound, lambda, prec);
  arf_ceil(ubound, ubound);
  res = arf_is_finite(ubound);

  if (res)
    {
      arf_get_fmpz_2exp(a, b,  ubound);
      fmpz_sub_si(b, b, prec/2 - THETA_DER_LOSS);
      mag_set_fmpz_2exp_fmpz(error, a, b);
    }

  arb_mat_clear(im);
  arb_clear(lambda);
  arb_clear(one);
  arf_clear(ubound);
  fmpz_clear(a);
  fmpz_clear(b);
  return res;
}

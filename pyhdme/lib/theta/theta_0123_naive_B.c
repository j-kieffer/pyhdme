
#include "theta.h"

int theta_0123_naive_B(fmpz_t B, const acb_mat_t tau, slong prec)
{
  arb_mat_t im;
  arb_t lambda, pi, q, b;
  arf_t sup;
  int res = 1;
  slong g = 2;

  arb_mat_init(im, g, g);
  arb_init(lambda);
  arb_init(pi);
  arb_init(q);
  arb_init(b);
  arf_init(sup);

  acb_mat_get_imag(im, tau);
  arb_mat_lambda(lambda, im, prec);

  if (!arb_is_positive(lambda)) res = 0;
  if (res)
    {
      arb_const_pi(pi, prec);
      arb_neg(q, pi);
      arb_mul(q, q, lambda, prec);
      arb_exp(q, q, prec);
      
      /* b = sqrt( (prec + 3 - 2log(1-q))/(pi lambda log_2(e)) ) - 1 */
      arb_sub_si(b, q, 1, prec);
      arb_neg(b, b);
      arb_log(b, b, prec);
      arb_mul_si(b, b, -2, prec);
      arb_add_ui(b, b, (ulong)prec + 3, prec);
      arb_div(b, b, pi, prec);
      arb_div(b, b, lambda, prec);
      /* Set pi to log_2(e) */
      arb_const_e(pi, prec);
      arb_log_base_ui(pi, pi, 2, prec);
      arb_div(b, b, pi, prec);
      
      arb_sqrt(b, b, prec);
      arb_sub_si(b, b, 1, prec);
      arb_get_ubound_arf(sup, b, prec);
      arf_ceil(sup, sup);
      arf_get_fmpz(B, sup, ARF_RND_NEAR);
    }
  
  arb_mat_clear(im);
  arb_clear(lambda);
  arb_clear(pi);
  arb_clear(q);
  arb_clear(b);
  arf_clear(sup);
  return res;
}

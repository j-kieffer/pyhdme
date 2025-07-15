
#include "theta.h"

int theta_use_naive(const acb_mat_t tau, slong prec)
{
  int res;
  arb_mat_t im;
  arb_t t;
  arb_t tol;

  arb_mat_init(im, 2, 2);
  arb_init(t);
  arb_init(tol);

  arb_one(tol);
  arb_mul_2exp_si(tol, tol, THETA_NEWTON_TOL_EXP);
  
  acb_mat_get_imag(im, tau);
  /* Add some tolerance to y1 */
  arb_one(t);
  arb_add(t, t, tol, prec);
  arb_mul(t, t, acb_imagref(acb_mat_entry(tau, 0, 0)), prec);
  
  res = siegel_is_weakly_reduced(tau, tol, prec);
  arb_mul_si(t, t, THETA_NEWTON_Y1, prec);
  arb_sub_si(t, t, prec, prec);
    res = res && arb_is_positive(t);
    
  arb_mat_clear(im);
  arb_clear(t);
  arb_clear(tol);
  
  return res;
}

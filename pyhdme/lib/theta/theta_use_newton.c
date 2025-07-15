
#include "theta.h"

int theta_use_newton(const acb_mat_t tau,  slong prec)
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
  arb_set(t, acb_imagref(acb_mat_entry(tau, 0, 0)));
  arb_sub_si(t, t, THETA_NEWTON_Y2MAX, prec);
  
  res = siegel_is_weakly_reduced(tau, tol, prec)
    && arb_is_negative(t);
  
  arb_mat_clear(im);
  arb_clear(t);
  arb_clear(tol);
  
  return res;
}


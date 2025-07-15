
#include "siegel.h"


int siegel_is_weakly_reduced(const acb_mat_t z, const arb_t tol, slong prec)
{
  slong g = acb_mat_nrows(z);
  arb_mat_t im;
  arb_t t;
  arb_t target;
  int res = 1;

  if (g != 2)
    {
      flint_printf("No notion of weak reduction for g = %wd\n", g);
      flint_abort();
    }

  arb_mat_init(im, g, g);
  arb_init(t);
  arb_init(target);

  acb_mat_get_imag(im, z);
  
  res = res && siegel_is_real_reduced(z, tol, prec)
    && arb_mat_is_minkowski_reduced(im, tol, prec); /* Actually stronger than needed */

  /* y1 \geq sqrt(3)/2, with tolerance */
  arb_add_si(t, tol, 1, prec);
  arb_mul(t, t, acb_imagref(acb_mat_entry(z, 0, 0)), prec);
  arb_set_si(target, 3);
  arb_sqrt(target, target, prec);
  arb_mul_2exp_si(target, target, -1);
  res = res && arb_ge(t, target);

  /* |z_1|, |z_2| \geq 1, with tolerance */
  arb_one(target);
  arb_sub(target, target, tol, prec);
  acb_abs(t, acb_mat_entry(z, 0, 0), prec);
  res = res && arb_ge(t, target);

  acb_abs(t, acb_mat_entry(z, 1, 1), prec);
  res = res && arb_ge(t, target);

  arb_mat_clear(im);
  arb_clear(t);
  arb_clear(target);

  return res;
}

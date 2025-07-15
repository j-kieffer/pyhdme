
#include "siegel.h"

void siegel_fundamental_domain_randtest(acb_mat_t z, flint_rand_t state, slong prec)
{
  arb_t tol;
  arb_t zero;
  acb_mat_t w;
  fmpz_mat_t m;
  slong g = 2;
  slong tol_bits = prec/5;
  int res = 0;

  arb_init(tol); /* Small tolerance during reduction */
  arb_init(zero);
  arb_one(tol);
  arb_mul_2exp_si(tol, tol, -tol_bits);
  acb_mat_init(w, g, g);
  fmpz_mat_init(m, 2*g, 2*g);

  while (!res)
    {
      siegel_halfspace_randtest(w, state, prec);
      res = siegel_fundamental_domain(w, m, w, tol, prec);
      if (res) res = siegel_is_in_fundamental_domain(w, zero, prec);
    }
  acb_mat_set(z, w);

  arb_clear(tol);
  arb_clear(zero);
  acb_mat_clear(w);
  fmpz_mat_clear(m);
}


#include "siegel.h"

void
arb_mat_randtest_sym_precise(arb_mat_t r, flint_rand_t state, slong prec, slong mag_bits)
{
  slong g = arb_mat_nrows(r);
  arb_mat_t rt;

  arb_mat_init(rt, g, g);

  arb_mat_randtest_precise(r, state, prec, mag_bits);
  arb_mat_transpose(rt, r);
  arb_mat_add(r, rt, r, prec);
  arb_mat_scalar_mul_2exp_si(r, r, -1);

  arb_mat_clear(rt);
}

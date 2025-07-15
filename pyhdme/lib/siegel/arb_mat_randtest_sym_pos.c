
#include "siegel.h"

void
arb_mat_randtest_sym_pos(arb_mat_t r, flint_rand_t state, slong prec)
{
  arb_mat_randtest_sym_precise(r, state, prec, 0);
  arb_mat_exp(r, r, prec);
}

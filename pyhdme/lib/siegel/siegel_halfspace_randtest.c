
#include "siegel.h"

void
siegel_halfspace_randtest(acb_mat_t z, flint_rand_t state, slong prec)
{
  arb_mat_t re, im;
  slong re_mag_bits = 0;
  
  slong g = acb_mat_ncols(z);
  prec = FLINT_MAX(20, prec);

  arb_mat_init(re, g, g);
  arb_mat_init(im, g, g);

  arb_mat_randtest_sym_precise(re, state, prec, re_mag_bits);
  arb_mat_randtest_sym_pos(im, state, prec);

  acb_mat_set_arb_arb(z, re, im);
  
  arb_mat_clear(re);
  arb_mat_clear(im);
}

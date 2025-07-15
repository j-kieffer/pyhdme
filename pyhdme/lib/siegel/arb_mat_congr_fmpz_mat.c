
#include "siegel.h"

void
arb_mat_congr_fmpz_mat(arb_mat_t r, const fmpz_mat_t u, const arb_mat_t m, slong prec)
{
  arb_mat_t x, u_arb;
  slong g = arb_mat_ncols(m);

  arb_mat_init(x, g, g);
  arb_mat_init(u_arb, g, g);
  arb_mat_set_fmpz_mat(u_arb, u);
  
  arb_mat_set(x, u_arb);
  arb_mat_transpose(x, x);
  arb_mat_mul(x, m, x, prec);
  arb_mat_mul(x, u_arb, x, prec);

  arb_mat_set(r, x);

  arb_mat_clear(x);
  arb_mat_clear(u_arb);
}
  
  
  

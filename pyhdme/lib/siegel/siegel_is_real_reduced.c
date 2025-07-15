
#include "siegel.h"

int
siegel_is_real_reduced(const acb_mat_t z, const arb_t tol, slong prec)
{
  arb_t bound;
  arb_t abs;
  int res = 1;
  slong g = acb_mat_nrows(z);
  int i, j;

  arb_init(bound);
  arb_init(abs);
  
  arb_one(bound);
  arb_mul_2exp_si(bound, bound, -1);
  arb_add(bound, bound, tol, prec);

  for (i = 0; i < g; i++)
    {
      for (j = 0; j < g; j++)
	{
	  arb_abs(abs, acb_realref(acb_mat_entry(z, i, j)));
	  res = res && arb_lt(abs, bound);
	}
    }

  arb_clear(bound);
  arb_clear(abs);
  return res;
}

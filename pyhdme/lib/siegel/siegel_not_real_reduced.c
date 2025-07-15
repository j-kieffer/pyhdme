
#include "siegel.h"

int
siegel_not_real_reduced(const acb_mat_t z, slong prec)
{
  arb_t abs;
  int res = 0;
  slong g = acb_mat_nrows(z);
  int i, j;

  arb_init(abs);

  for (i = 0; i < g; i++)
    {
      for (j = 0; j < g; j++)
	{
	  arb_abs(abs, acb_realref(acb_mat_entry(z, i, j)));
	  arb_mul_si(abs, abs, 2, prec);
	  arb_sub_si(abs, abs, 1, prec);
	  res = res || arb_is_positive(abs);
	}
    }
  
  arb_clear(abs);
  return res;
}

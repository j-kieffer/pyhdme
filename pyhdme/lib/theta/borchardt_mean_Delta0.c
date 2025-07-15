
#include "theta.h"

void
borchardt_mean_Delta0(arb_t Delta0, acb_srcptr a, slong prec)
{
  int i, j;
  acb_t diff;
  arf_t abs;
  arf_t Delta0_arf;

  acb_init(diff);
  arf_init(abs);
  arf_init(Delta0_arf);

  arf_zero(Delta0_arf);

  for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 4; j++)
	{
	  if (i != j)
	    {
	      acb_sub(diff, &a[i], &a[j], prec);
	      acb_get_abs_ubound_arf(abs, diff, prec);
	      arf_add(Delta0_arf, Delta0_arf, abs, prec, ARF_RND_CEIL);
	    }
	}
    }

  arb_set_arf(Delta0, Delta0_arf);

  acb_clear(diff);
  arf_clear(abs);
  arf_clear(Delta0_arf);
}

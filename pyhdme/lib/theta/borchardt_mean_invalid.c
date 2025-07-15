
#include "theta.h"

/* It is certain that the entries of the vector a do not fit in a half
   plane seen from the origin. */

/* The new version does not assume real 1 is part of the values: we
   detect if the excluded intervals cover S^1. Todo? We could also use
   this function to come up with the right angle to use in
   borchardt_mean_m0. */
int borchardt_mean_invalid(acb_srcptr a, slong prec)
{
  arf_struct ints[16];
  arf_t u, l;
  slong j, k;
  int res = 0;

  arf_init(u);
  arf_init(l);
  for (k = 0; k < 16; k++) arf_init(&ints[k]);
  
  for (k = 0; k < 4; k++) borchardt_excl_half_planes(&ints[4*k], &a[k], prec);
  /* Debug */
  /* flint_printf("(borchardt_mean_invalid) Values:\n");
  for (k = 0; k < 4; k++)
    {
      acb_printd(&a[k], 30); flint_printf("\n");
    }
  flint_printf("(borchardt_mean_invalid) Intervals:\n");
  for (k = 0; k < 16; k++)
    {
      arf_printd(&ints[k], 30); flint_printf("\n");
      } */
  
  arf_set_si(l, -1);
  for (k = 0; k < 4; k++) arf_max(l, l, &ints[4*k]);
  arf_set_si(u, 1);
  for (k = 0; k < 4; k++) arf_min(u, u, &ints[4*k+3]);

  for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
	{
	  if (arf_cmp(&ints[4*k+1], l) < 0)
	    {
	      arf_max(l, l, &ints[4*k+2]);
	    }
	}
    }

  if (arf_cmp(l, u) > 0) res = 1;

  for (k = 0; k < 16; k++) arf_clear(&ints[k]);
  arf_clear(u);
  arf_clear(l);
  return res;
}

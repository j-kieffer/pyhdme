
#include "theta.h"

int
borchardt_mean_quad_conv_is_reached(acb_srcptr a, slong prec)
{
  arb_t bound, absdiff;
  acb_t diff;
  int res = 1;
  int i;

  arb_init(bound);
  arb_init(absdiff);
  acb_init(diff);

  acb_abs(bound, &a[0], prec);
  arb_div_si(bound, bound, 7, prec);

  for (i = 1; i < 4; i++)
    {
      acb_sub(diff, &a[i], &a[0], prec);
      acb_abs(absdiff, diff, prec);
      res = res && arb_lt(absdiff, bound);
    }

  arb_clear(bound);
  arb_clear(absdiff);
  acb_clear(diff);
  
  return res;
}

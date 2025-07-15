
#include "igusa.h"

/* All three fundamental square theta quotients should have positive
   real part. */
int thomae_discard(acb_srcptr th2, slong prec)
{
  int res = 0;
  slong k;
  for (k = 1; k < 4; k++)
    {
      if (arb_is_nonpositive(acb_realref(&th2[k]))) res = 1;
    }
  res = res || theta2_invalid(th2, prec);
  return res;
}

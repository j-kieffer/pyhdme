
#include "theta.h"

/* If z is nonzero and intersects the negative real axis, multiply by
   -1 and output the root with positive imaginary part;
   otherwise, default root */

void borchardt_sqrt(acb_t r, const acb_t z, slong prec)
{
  acb_t res;
  acb_init(res);
  
  if (arb_contains_zero(acb_imagref(z))
      && arb_is_negative(acb_realref(z)))
    {
      acb_neg(res, z);
      acb_sqrt(res, res, prec);
      acb_mul_onei(res, res);
    }
  else acb_sqrt(res, z, prec);
  
  acb_set(r, res);
  acb_clear(res);
}

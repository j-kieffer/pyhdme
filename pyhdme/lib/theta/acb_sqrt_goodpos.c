
#include "theta.h"

int
acb_sqrt_goodpos(acb_t r, const acb_t z, slong prec)
{
  acb_sqrt(r, z, prec);
  return arb_is_positive(acb_realref(r));
}

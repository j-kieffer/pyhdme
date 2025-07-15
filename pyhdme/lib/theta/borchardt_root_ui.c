
#include "theta.h"

void borchardt_root_ui(acb_t r, const acb_t z, ulong e, slong prec)
{
  acb_t res;
  acb_t scal;
  acb_init(res);
  acb_init(scal);
  
  if (arb_contains_zero(acb_imagref(z))
      && arb_is_negative(acb_realref(z)))
    {
      acb_neg(res, z);
      acb_root_ui(res, res, e, prec);
      acb_set_si(scal, e);
      acb_inv(scal, scal, prec);
      acb_exp_pi_i(scal, scal, prec);
      acb_mul(res, res, scal, prec);
    }
  else acb_root_ui(res, z, e, prec);
  
  acb_set(r, res);
  acb_clear(res);
  acb_clear(scal);
}

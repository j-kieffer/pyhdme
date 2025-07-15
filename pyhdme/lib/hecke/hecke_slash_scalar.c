
#include "hecke.h"

void hecke_slash_scalar(acb_t im, const acb_t stardet, const acb_t val,
			slong k, slong prec)
{
  acb_t scal;

  acb_init(scal);
  
  acb_pow_si(scal, stardet, -k, prec);
  acb_mul(im, scal, val, prec);
  
  acb_clear(scal);
}

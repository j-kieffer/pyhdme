
#include "hecke.h"

void hecke_normalize_entry(hecke_t H, slong k,
			   fmpz* I, slong norm_ind, slong prec)
{  
  acb_t s;
  acb_t c;
  slong weights[4] = IGUSA_WEIGHTS;

  acb_init(s);
  acb_init(c);
  
  cov_find_rescaling(s, hecke_I_tau(H), I, 4, weights, prec);
  acb_inv(s, s, prec);
  acb_mul_si(s, s, norm_ind, prec);
  acb_mul_si(s, s, 6, prec);
  
  acb_inv(c, hecke_stardet(H, k), prec);
  acb_mul(c, c, s, prec);
  cov_rescale(hecke_I_norm(H, k), hecke_I(H, k), c, 4, weights, prec);
  
  acb_clear(s);
  acb_clear(c);
}

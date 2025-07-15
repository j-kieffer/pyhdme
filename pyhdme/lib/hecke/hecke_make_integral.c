
#include "hecke.h"

void hecke_make_integral(hecke_t H, fmpz* I, slong prec)
{
  slong k;
  acb_t s;
  acb_t c;
  slong weights[4] = IGUSA_WEIGHTS;

  acb_init(s);
  acb_init(c);
  
  cov_find_rescaling(s, hecke_I_tau(H), I, 4, weights, prec);
  acb_inv(s, s, prec);
  acb_mul_si(s, s, hecke_norm_ind(H), prec);
  acb_mul_si(s, s, 6, prec);

  for (k = 0; k < hecke_nb(H); k++)
    {
      acb_inv(c, hecke_stardet(H, k), prec);
      acb_mul(c, c, s, prec);
      cov_rescale(hecke_I_norm(H, k), hecke_I(H, k), c, 4, weights, prec);
    }
  
  acb_clear(s);
  acb_clear(c);
}

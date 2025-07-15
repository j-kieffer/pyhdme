
#include "hecke.h"

slong hecke_integral_highprec(const hecke_t H, slong prec)
{
  acb_t m;
  arf_t rad;
  slong res;
  slong j, k;

  acb_init(m);
  arf_init(rad);
  res = 0;

  for (j = 0; j < 4; j++)
    {
      acb_one(m);
      for (k = 0; k < hecke_nb(H); k++)
	{
	  acb_mul(m, m, &hecke_I_norm(H, k)[j], prec);
	  acb_mul_si(m, m, 2, prec);
	}
      acb_get_abs_ubound_arf(rad, m, prec);
      res = FLINT_MAX(res, arf_abs_bound_lt_2exp_si(rad));
    }
  acb_clear(m);
  arf_clear(rad);
  return 100*((prec + res)/90 + 1);
}

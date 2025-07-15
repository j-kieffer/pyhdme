
#include "hecke.h"

int hecke_has_integral_precision(hecke_t H, slong prec)
{
  arf_t rad;
  int res = 1;
  slong bits;
  slong k, j;
  
  arf_init(rad);

  for (k = 0; k < hecke_nb(H); k++)
    {
      for (j = 0; j < 4; j++)
	{
	  acb_get_rad_ubound_arf(rad, &hecke_I_norm(H, k)[j], prec);
	  bits = arf_abs_bound_lt_2exp_si(rad);
	  if (bits > -2)
	    {
	      res = 0;
	      break;
	    }
	}
      if (!res) break;
    }

  arf_clear(rad);
  return res;
}

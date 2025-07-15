
#include "hecke.h"

int hecke_set_I_fmpz_hilbert(hecke_t H, fmpz* I, slong delta, slong prec)
{
  int res;
  int v = get_hecke_verbose();

  res = hecke_set_I_fmpz(H, I, prec);
  if (res)
    {
      res = hilbert_inverse(hecke_t1t2(H), hecke_eta(H),
			    hecke_tau(H), delta, prec);
      if (v && !res)
	{
	  flint_printf("(hecke_set_I_fmpz_hilbert) Warning: Hilbert inversion failed\n");
	}
    }  
  return res;
}

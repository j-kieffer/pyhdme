
#include "igusa.h"

void cov_normalize_fmpz_wt1(fmpz* I, fmpz* S, slong nb)
{
  fmpz_t g;
  slong j;
  
  fmpz_init(g);
  _fmpz_vec_content(g, S, nb);
  _fmpz_vec_scalar_divexact_fmpz(I, S, nb, g);
  
  /* Check sign */
  for (j = 0; j < nb; j++)
    {
      if (!fmpz_is_zero(&I[j]))
	{
	  if (fmpz_cmp_si(&I[j], 0) < 0) _fmpz_vec_neg(I, I, nb); 
	  break;
	}
    }
  
  fmpz_clear(g);
}

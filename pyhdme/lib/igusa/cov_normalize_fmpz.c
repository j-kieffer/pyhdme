
#include "igusa.h"

void cov_normalize_fmpz(fmpz* I, fmpz* S, slong nb, slong* weights)
{
  fmpz_factor_t fac;
  slong k;

  fmpz_factor_init(fac);
  cov_factors(fac, S, nb);

  for (k = 0; k < nb; k++) fmpz_set(&I[k], &S[k]);
  
  for (k = 0; k < cov_factor_nb(fac); k++)
    {
      while (cov_divisible_fmpz(I, cov_factor_p(fac, k), nb, weights))
	{
	  cov_divexact_fmpz(I, I, cov_factor_p(fac, k), nb, weights);
	}
    }  
  
  /* Update a sign in case of odd weight */
  for (k = 0; k < nb; k++)
    {
      if (!fmpz_is_zero(&I[k]) && (weights[k] %2) == 1)
	{
	  if (fmpz_cmp_si(&I[k], 0) < 0) cov_rescale_fmpz_si(I, I, -1, nb, weights); 
	  break;
	}
    }

  fmpz_factor_clear(fac);
}

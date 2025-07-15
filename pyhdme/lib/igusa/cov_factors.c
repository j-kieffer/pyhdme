
#include "igusa.h"

static void cov_factors_additional(fmpz_factor_t fac, const fmpz_t g, fmpz* I, slong nb)
{
  fmpz_t r;
  slong k;
  int addg = 1;

  fmpz_init(r);

  for (k = 0; k < nb; k++)
    {
      if (!fmpz_is_zero(&I[k]))
	{
	  fmpz_remove(r, &I[k], g);
	  fmpz_gcd(r, g, r);
	  if (!fmpz_is_one(r))
	    {
	      cov_factors_additional(fac, r, I, nb);
	      fmpz_divexact(r, g, r);
	      cov_factors_additional(fac, r, I, nb);
	      addg = 0;
	      break;
	    }
	}
    }

  if (addg) _fmpz_factor_append(fac, g, 1);
  fmpz_clear(r);
}

void cov_factors(fmpz_factor_t fac, fmpz* I, slong nb)
{
  fmpz_t g;
  slong k;
  int stop;

  fmpz_init(g);
  fmpz_zero(g);
  
  for (k = 0; k < nb; k++)
    {
      fmpz_gcd(g, g, &I[k]);
    }
  if (fmpz_is_zero(g))
    {
      flint_printf("(cov_factors) Error: all entries are zero\n");
      fflush(stdout);
      flint_abort();
    }

  stop = fmpz_factor_smooth(fac, g, COV_FACTOR_BITS, 1);
  for (k = 0; k < cov_factor_nb(fac) + (stop ? 0 : -1); k++) fmpz_remove(g, g, cov_factor_p(fac, k));
  if (!stop) cov_factors_additional(fac, g, I, nb);

  fmpz_clear(g);
}

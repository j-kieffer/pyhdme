
#include "polynomials.h"

void pol_remove_factor_Q(fmpz_poly_t r, const fmpz_poly_t pol,
			 const fmpz_poly_t fac, slong mult)
{
  
  fmpz_poly_t rem;
  slong k;

  fmpz_poly_init(rem);
  
  fmpz_poly_set(r, pol);
  for (k = 0; k < mult; k++)
    {
      fmpz_poly_divrem(r, rem, r, fac);
      if (!fmpz_poly_is_zero(rem))
	{
	  flint_printf("(pol_remove_factor_Q) Not divisible\n");
	  fflush(stdout);
	  flint_abort();
	}      
    }

  fmpz_poly_clear(rem);
}


#include "polynomials.h"

/* See also pol_remove_root_Q */

void pol_remove_root_Fp(fmpz_mod_poly_t r, const fmpz_mod_poly_t pol,
			const fmpz_t root, slong mult, const fmpz_mod_ctx_t ctx)
{  
  fmpz_mod_poly_t fac;
  fmpz_mod_poly_t rem;
  fmpz_t c;
  slong k;

  fmpz_mod_poly_init(fac, ctx);
  fmpz_mod_poly_init(rem, ctx);
  fmpz_init(c);

  /* Set fac to a degree 1 polynomial with given root */
  fmpz_mod_poly_set_coeff_ui(fac, 1, 1, ctx);
  fmpz_neg(c, root);
  fmpz_mod_poly_set_coeff_fmpz(fac, 0, c, ctx);
  
  fmpz_neg(fmpz_poly_get_coeff_ptr(fac, 0), fmpz_poly_get_coeff_ptr(fac, 0));  

  fmpz_mod_poly_set(r, pol, ctx);
  for (k = 0; k < mult; k++)
    {
      fmpz_mod_poly_divrem(r, rem, r, fac, ctx);
      if (!fmpz_poly_is_zero(rem))
	{
	  flint_printf("(pol_remove_root_Fp) Not divisible\n");
	  fflush(stdout);
	  flint_abort();
	}      
    }

  fmpz_mod_poly_clear(fac, ctx);
  fmpz_mod_poly_clear(rem, ctx);
  fmpz_clear(c);
}

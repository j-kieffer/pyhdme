
#include "polynomials.h"

void pol_remove_root_Q(fmpz_poly_t r, const fmpz_poly_t pol,
		       const fmpq_t root, slong mult)
{
  fmpz_poly_t fac;

  fmpz_poly_init(fac);

  /* Set fac to a degree 1 polynomial with given root */
  fmpz_poly_set_coeff_fmpz(fac, 1, fmpq_denref(root));
  fmpz_poly_set_coeff_fmpz(fac, 0, fmpq_numref(root));
  fmpz_neg(fmpz_poly_get_coeff_ptr(fac, 0), fmpz_poly_get_coeff_ptr(fac, 0));

  pol_remove_factor_Q(r, pol, fac, mult);

  fmpz_poly_clear(fac);
}

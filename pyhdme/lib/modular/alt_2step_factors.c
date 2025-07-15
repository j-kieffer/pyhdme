
#include "modular.h"

void alt_2step_factors(slong* nb, fmpz_poly_struct* factors, slong* mults,
		       const modeq_t E, slong ell)
{
  slong nb_factors = 0;
  slong max_nb_factors = fmpz_poly_degree(modeq_equation(E));
  fmpz_poly_struct* all_factors;
  slong* exps;
  fmpz_poly_t fac;
  slong d;
  slong k;

  all_factors = flint_malloc(max_nb_factors * sizeof(fmpz_poly_struct));
  for (k = 0; k < max_nb_factors; k++) fmpz_poly_init(&all_factors[k]);
  exps = flint_malloc(max_nb_factors * sizeof(slong));
  fmpz_poly_init(fac);
  
  pol_factor_Q(&nb_factors, all_factors, exps, modeq_equation(E));
  *nb = 0;
  
  /* Isolate all factors of degree not 1 and dividing l+1 */
  for (k = 0; k < nb_factors; k++)
    {
      d = fmpz_poly_degree(&all_factors[k]);
      if (d > 1 && (ell+1) % d == 0)
	{
	  fmpz_poly_set(fac, &all_factors[k]);
	  fmpz_poly_primitive_part(fac, fac);
	  fmpz_poly_set(&factors[*nb], fac);
	  mults[*nb] = exps[k];
	  (*nb)++;
	}
    }

  for (k = 0; k < max_nb_factors; k++) fmpz_poly_clear(&all_factors[k]);
  flint_free(all_factors);
  flint_free(exps);
  fmpz_poly_clear(fac);
}

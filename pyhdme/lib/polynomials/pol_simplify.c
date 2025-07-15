
#include "polynomials.h"

void pol_simplify(fmpz_poly_struct* num_vec, fmpz_t den, slong degree, slong nb)
{
  fmpz_t gcd;
  fmpz_t coeff;
  slong gcd_bits;
  slong i, k;
  int v = get_pol_verbose();
  
  fmpz_init(gcd);
  fmpz_init(coeff);
  
  fmpz_set(gcd, den);
  for (i = 0; i < nb; i++)
    {
      for (k = 0; k <= degree ; k++)
	{
	  fmpz_poly_get_coeff_fmpz(coeff, &num_vec[i], k);
	  fmpz_gcd(gcd, gcd, coeff);
	}
    }
  
  fmpz_divexact(den, den, gcd);
  for (i = 0; i < nb; i++)
    {
      fmpz_poly_scalar_divexact_fmpz(&num_vec[i], &num_vec[i], gcd);
    }
  
  gcd_bits = fmpz_bits(gcd);
  if (v)
    {
      flint_printf("(pol_simplify) Simplify result by %wd-bit integer\n", gcd_bits);
      fmpz_set_si(coeff, 2);
      flint_printf("(pol_simplify) Valuation at 2: %wd\n", fmpz_remove(gcd ,gcd, coeff));
      fmpz_set_si(coeff, 3);
      flint_printf("(pol_simplify) Valuation at 3: %wd\n", fmpz_remove(gcd, gcd, coeff));
      fmpz_set_si(coeff, 5);
      flint_printf("(pol_simplify) Valuation at 5: %wd\n", fmpz_remove(gcd, gcd, coeff));
      fmpz_set_si(coeff, 7);
      flint_printf("(pol_simplify) Valuation at 7: %wd\n", fmpz_remove(gcd, gcd, coeff));
      flint_printf("(pol_simplify) Rest: %wd bits\n", fmpz_bits(gcd));
    }    
  
  fmpz_clear(gcd);
  fmpz_clear(coeff);
}

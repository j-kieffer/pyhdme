
#include "hilbert.h"

void hilbert_conjugate(fmpz_poly_t xbar, const fmpz_poly_t x, slong delta)
{
  fmpz_poly_t z;
  fmpz_t c, temp;
  
  fmpz_poly_init(z);
  fmpz_init(c);
  fmpz_init(temp);

  fmpz_poly_get_coeff_fmpz(c, x, 0);
  fmpz_poly_get_coeff_fmpz(temp, x, 1);
  if (delta % 2 == 1)
    {
      fmpz_add(c, c, temp);
    }			
  fmpz_poly_set_coeff_fmpz(z, 0, c);
  fmpz_neg(c, temp);
  fmpz_poly_set_coeff_fmpz(z, 1, c);

  fmpz_poly_set(xbar, z);
  
  fmpz_poly_clear(z);
  fmpz_clear(c);
  fmpz_clear(temp);  
}

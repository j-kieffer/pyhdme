
#include "hilbert.h"

void hilbert_transform_randtest(fmpz_poly_mat_t m, flint_rand_t state, slong bits)
{
  fmpz_mat_t x;
  fmpz_t det;
  slong j, k;
  
  fmpz_mat_init(x, 2, 2);
  fmpz_init(det);

  fmpz_one(det);
  fmpz_mat_randdet(x, state, det);
  fmpz_mat_randops(x, state, 10);
  
  fmpz_poly_mat_zero(m);
  for (j = 0; j < 2; j++)
    {
      for (k = 0; k < 2; k++)
	{
	    fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(m, j, k),
				     0, fmpz_mat_entry(x, j, k));
	}
    }

  fmpz_mat_clear(x);
  fmpz_clear(det);
}


#include "siegel.h"

void fmpz_mat_get_b(fmpz_mat_t b, const fmpz_mat_t m)
{
  slong g = fmpz_mat_half_dim(m);
  slong j, k;
  for (j = 0; j < g; j++)
    {
      for (k = 0; k < g; k++)
	{
	  fmpz_set(fmpz_mat_entry(b, j, k),
		   fmpz_mat_entry(m, j, k+g));
	}
    }
}

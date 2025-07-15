
#include "siegel.h"

void
arb_mat_randtest_precise(arb_mat_t r, flint_rand_t state, slong prec, slong mag_bits)
{
  slong nrows = arb_mat_nrows(r);
  slong ncols = arb_mat_ncols(r);
  slong i, j;
  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < ncols; j++)
	{
	  arb_randtest_precise(arb_mat_entry(r, i, j), state, prec, mag_bits);
	}
    }
}

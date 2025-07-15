
#include "siegel.h"

void acb_mat_set_window(acb_mat_t z, slong j, slong k, const acb_mat_t w)
{
  slong u, v;

  if (j + acb_mat_nrows(w) > acb_mat_nrows(z)
      || k + acb_mat_ncols(w) > acb_mat_ncols(z))
    {
      flint_printf("(acb_mat_set_windows) Wrong dimensions\n");
      fflush(stdout);
      flint_abort();
    }
  
  for (u = 0; u < acb_mat_nrows(w); u++)
    {
      for (v = 0; v < acb_mat_ncols(w); v++)
	{
	  acb_set(acb_mat_entry(z, j+u, k+v),
		  acb_mat_entry(w, u, v));
	}
    }
}

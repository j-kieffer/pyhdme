
#include "siegel.h"

int
siegel_not_in_fundamental_domain(const acb_mat_t z, slong prec)
{
  slong g = acb_mat_nrows(z);
  int res = 0;
  int j;

  fmpz_mat_t test;
  acb_mat_t star;
  arb_mat_t im;
  acb_t det;
  arb_t absdet;

  fmpz_mat_init(test, 2*g, 2*g);
  acb_mat_init(star, g, g);
  arb_mat_init(im, g, g);
  acb_init(det);
  arb_init(absdet);

  acb_mat_get_imag(im, z);
  res = siegel_not_real_reduced(z, prec)
    || arb_mat_not_minkowski_reduced(im, prec);

  if (!res)
    {
      /* Test matrices */
      for (j = 0; j < siegel_nb_test_matrices(g); j++)
	{
	  siegel_test_matrix(test, j);
	  siegel_star(star, test, z, prec);
	  acb_mat_det(det, star, prec);
	  acb_abs(absdet, det, prec);
	  arb_sub_si(absdet, absdet, 1, prec);
	  if (arb_is_negative(absdet))
	    {
	      res = 1;
	      break; /* for loop */
	    }
	}
    }

  fmpz_mat_clear(test);
  acb_mat_clear(star);
  arb_mat_clear(im);
  acb_clear(det);
  arb_clear(absdet);
  return res;
}


#include "siegel.h"


int
siegel_is_in_fundamental_domain(const acb_mat_t z, const arb_t tol, slong prec)
{
  slong g = acb_mat_nrows(z);
  int res = 1;
  int j;

  fmpz_mat_t test;
  acb_mat_t star;
  arb_mat_t im;
  acb_t det;
  arb_t absdet;
  arb_t one_minus_eps;

  fmpz_mat_init(test, 2*g, 2*g);
  acb_mat_init(star, g, g);
  arb_mat_init(im, g, g);
  acb_init(det);
  arb_init(absdet);
  arb_init(one_minus_eps);

  arb_set_si(one_minus_eps, 1);
  arb_sub(one_minus_eps, one_minus_eps, tol, prec);

  acb_mat_get_imag(im, z);
  res = siegel_is_real_reduced(z, tol, prec)
    &&  arb_mat_is_minkowski_reduced(im, tol, prec);
  
  if (res)
    {
      /* Test matrices */
      for (j = 0; j < siegel_nb_test_matrices(g); j++)
	{
	  siegel_test_matrix(test, j);
	  siegel_star(star, test, z, prec);
	  acb_mat_det(det, star, prec);
	  acb_abs(absdet, det, prec);
	  if (!arb_gt(absdet, one_minus_eps))
	    {
	      res = 0;
	      break; /* for loop */
	    }
	}
    }
  
  fmpz_mat_clear(test);
  acb_mat_clear(star);
  arb_mat_clear(im);
  acb_clear(det);
  arb_clear(absdet);
  arb_clear(one_minus_eps);
  return res;
}

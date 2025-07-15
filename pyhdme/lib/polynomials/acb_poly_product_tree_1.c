
#include "polynomials.h"

void acb_poly_product_tree_1(acb_poly_t P, acb_srcptr xi,
			     acb_srcptr yi, slong d, slong prec)
{
  acb_poly_struct* level;
  slong k, m, mmax;

  if (d <= 0)
    {
      acb_poly_one(P);
      return;
    }
  
  mmax = 1;
  while (mmax < d) mmax *= 2;
  level = flint_malloc(mmax * sizeof(acb_poly_struct));

  if (level == NULL)
    {
      flint_printf("(product_tree_1) Memory allocation fails\n");
      fflush(stdout);
      flint_abort();
    }
  for (k = 0; k < mmax; k++) acb_poly_init(&level[k]);

  m = 1;
  for (k = 0; k < d; k++)
    {
      acb_poly_set_coeff_acb(&level[k], 1, &xi[k]);
      acb_poly_set_coeff_acb(&level[k], 0, &yi[k]);
    }
  for (k = d; k < mmax; k++)
    {
      acb_poly_one(&level[k]);
    }
  while (m < mmax)
    {
      m *= 2;
      for (k = 0; k < mmax/m; k++)
	{
	  acb_poly_mul(&level[k], &level[2*k], &level[2*k+1], prec);
	}
    }
  acb_poly_set(P, &level[0]);

  for (k = 0; k < mmax; k++) acb_poly_clear(&level[k]);
  flint_free(level);
}


#include "polynomials.h"

void acb_poly_product_tree_2(acb_poly_t Q, acb_srcptr xi, acb_srcptr yi,
			     acb_srcptr zi, slong d, slong prec)
{
  acb_poly_struct* levelP;
  acb_poly_struct* levelQ;
  acb_poly_t prod1, prod2;
  slong k, m, mmax;

  if (d <= 0)
    {
      acb_poly_zero(Q);
      return;
    }
  
  mmax = 1;
  while (mmax < d) mmax *= 2;
  levelP = flint_malloc(mmax * sizeof(acb_poly_struct));
  levelQ = flint_malloc(mmax * sizeof(acb_poly_struct));

  if (levelP == NULL || levelQ == NULL)
    {
      flint_printf("(product_tree_1) Memory allocation fails\n");
      fflush(stdout);
      flint_abort();
    }

  acb_poly_init(prod1);
  acb_poly_init(prod2);
  for (k = 0; k < mmax; k++) acb_poly_init(&levelP[k]);
  for (k = 0; k < mmax; k++) acb_poly_init(&levelQ[k]);

  /* Set up level zero */
  m = 1;
  for (k = 0; k < d; k++)
    {
      acb_poly_set_coeff_acb(&levelP[k], 1, &xi[k]);
      acb_poly_set_coeff_acb(&levelP[k], 0, &yi[k]);
    }
  for (k = d; k < mmax; k++) acb_poly_one(&levelP[k]);
  for (k = 0; k < d; k++) acb_poly_set_acb(&levelQ[k], &zi[k]);
  for (k = d; k < mmax; k++) acb_poly_zero(&levelQ[k]);
  /* Next level from the previous one */
  while (m < mmax)
    {
      m *= 2;
      for (k = 0; k < mmax/m; k++)
	{
	  acb_poly_mul(prod1, &levelQ[2*k], &levelP[2*k+1], prec);
	  acb_poly_mul(prod2, &levelQ[2*k+1], &levelP[2*k], prec);
	  acb_poly_add(&levelQ[k], prod1, prod2, prec);
	  acb_poly_mul(&levelP[k], &levelP[2*k], &levelP[2*k+1], prec);
	}
    }
  acb_poly_set(Q, &levelQ[0]);
  
  acb_poly_clear(prod1);
  acb_poly_clear(prod2);
  for (k = 0; k < mmax; k++) acb_poly_clear(&levelP[k]);
  for (k = 0; k < mmax; k++) acb_poly_clear(&levelQ[k]);
  flint_free(levelP);
  flint_free(levelQ);
}



#include "hilbert.h"

/* Here I is the Streng covariants */

static void complete_from_G2(acb_t F6, acb_t F10, acb_t x, acb_t y,
			     const acb_t G2, acb_srcptr I, slong prec)
{
  acb_t temp;
  acb_init(temp);
  
  acb_div_si(F6, &I[1], 4, prec);
  acb_pow_si(temp, G2, 3, prec);
  acb_sub(F6, F6, temp, prec);
  acb_div_si(F6, F6, -864, prec);

  acb_div_si(F10, &I[2], n_pow(2,12), prec);

  acb_set(x, &I[3]);
  acb_div_si(x, x, n_pow(2,15), prec);
  acb_mul(y, F6, F6, prec);
  acb_mul_si(y, y, 3, prec);
  acb_mul(temp, G2, F10, prec);
  acb_mul_si(temp, temp, 2, prec);
  acb_sub(y, y, temp, prec);

  acb_clear(temp);
}

void gundlach_from_igusa(acb_ptr G, acb_srcptr I, slong delta, slong prec)
{
  acb_ptr S;
  acb_t G2, F6, F10;
  acb_t x, y;
  
  if (delta != 5)
    {
      flint_printf("(gundlach_from_igusa) Error: Gundlach invariants only implemented for discriminant 5\n");
      fflush(stdout);
      flint_abort();
    }

  acb_init(G2);
  acb_init(F6);
  acb_init(F10);
  acb_init(x);
  acb_init(y);
  S = _acb_vec_init(4);

  igusa_streng(S, I, prec);

  acb_div_si(G2, &S[0], 4, prec);
  borchardt_sqrt(G2, G2, prec); /* Possible sign error */
  
  /* Check sign by computing chi12 in two ways */
  complete_from_G2(F6, F10, x, y, G2, S, prec);
  if (!acb_overlaps(x, y))
    {
      acb_neg(G2, G2);
      complete_from_G2(F6, F10, x, y, G2, S, prec);
      if (!acb_overlaps(x, y))
	{
	  flint_printf("(gundlach_from_igusa) Could not compute corresponding Gundlach covariants\n");
	  fflush(stdout);
	  flint_abort();
	}
    }

  acb_set(&G[0], G2);
  acb_set(&G[1], F6);
  acb_set(&G[2], F10);
  
  acb_clear(G2);
  acb_clear(F6);
  acb_clear(F10);
  acb_clear(x);
  acb_clear(y);
  _acb_vec_clear(S, 4);
}

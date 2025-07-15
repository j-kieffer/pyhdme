
#include "igusa.h"

static void nonzero_indices(slong* inds, slong* k, fmpz* I, slong nb)
{
  slong j;
  *k = 0;
  for (j = 0; j < nb; j++)
    {
      if (!fmpz_is_zero(&I[j]))
	{
	  inds[*k] = j;
	  (*k)++;
	}
    }
}

static void extract_inds(slong* extr, slong* inds, slong k, slong* vec)
{
  slong j;
  for (j = 0; j < k; j++)
    {
      extr[j] = vec[inds[j]];
      if (vec[inds[j]] < 0)
	{
	  flint_printf("(cov_min_weight_combination) Error: negative weight\n");
	  fflush(stdout);
	  flint_abort();
	}
    }
}

static void xgcd_vec_si(slong* gcd, slong* coefs, slong* a, slong nb)
{
  slong j;
  ulong u, v;
  slong i;
  
  if (nb == 0)
    {
      flint_printf("(cov_min_weight_combination: xgcd_vec_si) nb = 0\n");
      fflush(stdout);
      flint_abort();
    }

  *gcd = a[0];
  coefs[0] = 1;
  /* gcd contains result for indices 0,...,j-1 */
  for (j = 1; j < nb; j++)
    {
      if (a[j] % (*gcd) == 0)
	{
	  coefs[j] = 0;
	}
      else if (*gcd >= a[j])
	{
	  *gcd = n_xgcd(&u, &v, *gcd, a[j]);
	  coefs[j] = -v;
	  for (i = 0; i < j; i++) coefs[i] *= u;
	}
      else
	{
	  *gcd = n_xgcd(&u, &v, a[j], *gcd);
	  coefs[j] = u;
	  for (i = 0; i < j; i++) coefs[i] *= (-v);
	}
    }
}

void cov_min_weight_combination(slong* wt, slong* exponents, fmpz* I,
				slong nb, slong* weights)
{
  slong* inds;
  slong* extr;
  slong* coefs;
  slong g;
  slong k;
  slong i;

  inds = flint_malloc(nb * sizeof(slong));
  extr = flint_malloc(nb * sizeof(slong));
  coefs = flint_malloc(nb * sizeof(slong));

  nonzero_indices(inds, &k, I, nb);
  if (k == 0)
    {
      flint_printf("(cov_min_weight_combination) Error: all entries are zero\n");
      fflush(stdout);
      flint_abort();
    }
  extract_inds(extr, inds, k, weights);
  xgcd_vec_si(&g, coefs, extr, k);
  
  *wt = g;
  for (i = 0; i < nb; i++) exponents[i] = 0;
  for (i = 0; i < k; i++) exponents[inds[i]] = coefs[i];

  flint_free(inds);
  flint_free(extr);
  flint_free(coefs);  
}

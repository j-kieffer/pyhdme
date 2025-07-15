
#include "modular.h"

int alt_2step_select_isog(slong* indices, const fmpz_poly_t factor, slong mult,
			  const hecke_t H, const modeq_ctx_t ctx)
{
  acb_t num, den;
  acb_poly_t pol;
  slong* seen;
  slong prec = hecke_prec(H);
  
  slong k, j;
  slong nb;
  slong d = fmpz_poly_degree(factor);
  int new;
  int res = 1;

  acb_init(num);
  acb_init(den);
  acb_poly_init(pol);
  seen = flint_malloc(hecke_nb(H) * sizeof(slong));
  
  acb_poly_set_fmpz_poly(pol, factor, prec);
  nb = 0;
  for (k = 0; k < hecke_nb(H); k++)
    {
      new = 1;
      /* If one of the ctx pairs, add 1 to "seen" and continue */
      for (j = 0; j < k; j++)
	{
	  if (modeq_ctx_is_pair(j, k, ctx))
	    {
	      new = 0;
	      seen[j]++;
	      break;
	    }
	}
      if (!new) continue;
      
      /* Compute associated root; is it a root of factor? */
      cov_mpoly_eval(den, modeq_ctx_den(ctx), hecke_I(H, k), modeq_ctx_ctx(ctx), prec);
      cov_mpoly_eval(num, modeq_ctx_num(ctx), hecke_I(H, k), modeq_ctx_ctx(ctx), prec);
      acb_div(num, num, den, prec);
      acb_poly_evaluate(num, pol, num, prec);

      if (acb_contains_zero(num)) /* Probable root */
	{
	  if (nb >= d) /* Too many roots, abort with res = 0 */
	    {
	      res = 0;
	      break;
	    }
	  indices[nb] = k;
	  seen[k] = 1;
	  nb++;
	}
    }

  if (res) /* We collected some roots, not too many. */
    {
      if (nb != d) res = 0;
      for (k = 0; k < nb; k++)
	{
	  if (seen[indices[k]] != mult) res = 0;
	}
    }
  
  acb_clear(num);
  acb_clear(den);
  acb_poly_clear(pol);
  flint_free(seen);
  return res;
}

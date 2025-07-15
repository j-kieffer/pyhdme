
#include "igusa.h"

void cov_eval_all_monomials(acb_ptr ev, acb_srcptr I, slong wt,
			    slong nb, slong* weights, slong prec)
{
  slong m = cov_nb_monomials(wt, nb, weights);
  slong* exps;
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t mon;
  slong k;

  exps = flint_malloc(m * nb * sizeof(slong));
  fmpz_mpoly_ctx_init(ctx, nb, ORD_LEX);
  fmpz_mpoly_init(mon, ctx);
  
  cov_all_exps(exps, wt, nb, weights);
  for (k = 0; k < m; k++)
    {
      cov_monomial(mon, &exps[nb * k], ctx);
      cov_mpoly_eval(&ev[k], mon, I, ctx, prec);
    }

  flint_free(exps);
  fmpz_mpoly_clear(mon, ctx);
  fmpz_mpoly_ctx_clear(ctx);
}

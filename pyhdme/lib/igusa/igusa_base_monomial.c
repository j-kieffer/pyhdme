
#include "igusa.h"

void igusa_base_monomial(fmpz_mpoly_t mon, slong wt, slong k,
			 const fmpz_mpoly_ctx_t ctx)
{
  slong exps[4];
  igusa_base_exps(exps, wt, k);
  cov_monomial(mon, exps, ctx);
}

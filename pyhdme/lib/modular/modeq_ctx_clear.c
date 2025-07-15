
#include "modular.h"

void modeq_ctx_clear(modeq_ctx_t ctx)
{
  slong j;
  slong m = MODEQ_MAX_NB_MONOMIALS;
  
  for (j = 0; j < m; j++)
    {
      fmpz_mpoly_clear(modeq_ctx_monomial(ctx, j), modeq_ctx_ctx(ctx));
    }
  fmpz_mpoly_clear(modeq_ctx_den(ctx), modeq_ctx_ctx(ctx));
  fmpz_mpoly_clear(modeq_ctx_num(ctx), modeq_ctx_ctx(ctx));
  
  flint_free(ctx->monomials);
  fmpz_mpoly_ctx_clear(modeq_ctx_ctx(ctx));
  flint_free(ctx->pairs);
}

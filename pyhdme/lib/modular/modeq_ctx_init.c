
#include "modular.h"

void modeq_ctx_init(modeq_ctx_t ctx)
{
  slong j;
  slong m = MODEQ_MAX_NB_MONOMIALS;
  
  fmpz_mpoly_ctx_init(modeq_ctx_ctx(ctx), 4, ORD_LEX);
  ctx->monomials = flint_malloc(m * sizeof(fmpz_mpoly_struct));
  for (j = 0; j < m; j++)
    {
      fmpz_mpoly_init(modeq_ctx_monomial(ctx, j), modeq_ctx_ctx(ctx));
    }
  fmpz_mpoly_init(modeq_ctx_den(ctx), modeq_ctx_ctx(ctx));
  fmpz_mpoly_init(modeq_ctx_num(ctx), modeq_ctx_ctx(ctx));

  ctx->pairs = flint_malloc(MODEQ_CTX_ALLOC * 2 * sizeof(slong));
  ctx->alloc_pairs = MODEQ_CTX_ALLOC;
  ctx->nb_pairs = 0;						     
}

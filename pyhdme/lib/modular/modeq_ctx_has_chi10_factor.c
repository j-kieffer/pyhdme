
#include "modular.h"

int modeq_ctx_has_chi10_factor(const modeq_ctx_t ctx)
{
  fmpz_mpoly_t chi10, r;
  slong exps[4] = {0,0,1,0};
  int res;

  fmpz_mpoly_init(chi10, modeq_ctx_ctx(ctx));
  fmpz_mpoly_init(r, modeq_ctx_ctx(ctx));  

  cov_monomial(chi10, exps, modeq_ctx_ctx(ctx));
  res = fmpz_mpoly_divides(r, modeq_ctx_num(ctx), chi10, modeq_ctx_ctx(ctx))
    && fmpz_mpoly_divides(r, modeq_ctx_den(ctx), chi10, modeq_ctx_ctx(ctx));
    
  fmpz_mpoly_clear(chi10, modeq_ctx_ctx(ctx));
  fmpz_mpoly_clear(r, modeq_ctx_ctx(ctx));
  return res;
}

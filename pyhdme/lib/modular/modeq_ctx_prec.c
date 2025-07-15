
#include "modular.h"

slong modeq_ctx_prec(slong prec)
{
  slong res;
  res = FLINT_MAX(MODEQ_CTX_MIN_PREC, n_sqrt(prec));
  return res;
}

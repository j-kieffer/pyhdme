
#include "hdme_data.h"

void hdme_data_evaluate_fmpq(fmpq_t ev, const fmpq_mpoly_t pol,
			     fmpq* vals, const fmpq_mpoly_ctx_t ctx)
{
  slong nb = fmpq_mpoly_ctx_nvars(ctx);
  fmpq** ptr;
  slong k;

  ptr = flint_malloc(nb * sizeof(fmpq*));
  for (k = 0; k < nb; k++) ptr[k] = &vals[k];

  fmpq_mpoly_evaluate_all_fmpq(ev, pol, (fmpq* const*) ptr, ctx);

  flint_free(ptr);
}

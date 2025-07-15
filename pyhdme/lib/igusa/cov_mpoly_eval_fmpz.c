
#include "igusa.h"

/* See also hdme_data_evaluate_fmpq */

void cov_mpoly_eval_fmpz(fmpz_t ev, const fmpz_mpoly_t pol, fmpz* I,
			 const fmpz_mpoly_ctx_t ctx)
{
  slong nb = fmpz_mpoly_ctx_nvars(ctx);
  fmpz** ptr;
  slong k;

  ptr = flint_malloc(nb * sizeof(fmpz*));
  for (k = 0; k < nb; k++) ptr[k] = &I[k];

  fmpz_mpoly_evaluate_all_fmpz(ev, pol, (fmpz* const*) ptr, ctx);

  flint_free(ptr);
}

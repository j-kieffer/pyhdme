
#include "igusa.h"

/* Fix minimum weight combination to one. Possible divisions by zero. */

void cov_normalize(acb_ptr I, acb_srcptr S, slong nb, slong* weights, slong prec)
{
  fmpz* fake;
  slong wt;
  slong* exps;
  acb_t r, c;
  slong j;

  fake = _fmpz_vec_init(nb);
  exps = flint_malloc(nb * sizeof(slong));
  acb_init(c);
  acb_init(r);

  for (j = 0; j < nb; j++) fmpz_one(&fake[j]);
  cov_min_weight_combination(&wt, exps, fake, nb, weights);

  acb_one(c);
  for (j = 0; j < nb; j++)
    {
      acb_pow_si(r, &S[j], exps[j], prec);
      acb_mul(c, c, r, prec);
    }
  acb_inv(c, c, prec);
  borchardt_root_ui(c, c, wt, prec);
  cov_rescale(I, S, c, nb, weights, prec);

  _fmpz_vec_clear(fake, nb);
  flint_free(exps);
  acb_clear(c);
  acb_clear(r);
}

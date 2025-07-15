
#include "igusa.h"

/* See also hdme_data_evaluate_acb */

void cov_mpoly_eval(acb_t ev, const fmpz_mpoly_t pol, acb_srcptr I,
		    const fmpz_mpoly_ctx_t ctx, slong prec)
{
  slong n = fmpz_mpoly_ctx_nvars(ctx);
  slong L = fmpz_mpoly_length(pol, ctx);
  slong* degrees = flint_malloc(n * sizeof(slong));
  slong j, k;
  acb_ptr* powers = flint_malloc(n * sizeof(acb_ptr));
  acb_t res, temp;
  fmpz_t coeff;
  slong exp;

  fmpz_mpoly_degrees_si(degrees, pol, ctx);
  for (k = 0; k < n; k++)
    {
      powers[k] = _acb_vec_init(degrees[k]+2);
      acb_one(&(powers[k][0]));
      for (j = 1; j <= degrees[k]; j++)
	{
	  acb_mul(&(powers[k][j]), &(powers[k][j-1]), &I[k], prec);
	}
    }
  acb_init(res);
  acb_init(temp);
  fmpz_init(coeff);

  acb_zero(res);
  for (j = 0; j < L; j++)
    {
      fmpz_mpoly_get_term_coeff_fmpz(coeff, pol, j, ctx);
      acb_set_fmpz(temp, coeff);
      for (k = 0; k < n; k++)
	{
	  exp = fmpz_mpoly_get_term_var_exp_si(pol, j, k, ctx);
	  acb_mul(temp, temp, &(powers[k][exp]), prec);
	}
      acb_add(res, res, temp, prec);
    }

  acb_set(ev, res);

  acb_clear(res);
  acb_clear(temp);
  fmpz_clear(coeff);
  for (k = 0; k < n; k++)
    {
      _acb_vec_clear(powers[k], degrees[k]+2);
    }
  flint_free(degrees);
  flint_free(powers);
}

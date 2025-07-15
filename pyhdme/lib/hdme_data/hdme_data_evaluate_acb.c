
#include "hdme_data.h"

void hdme_data_evaluate_acb(acb_t ev, const fmpq_mpoly_t pol, acb_srcptr vals,
			    const fmpq_mpoly_ctx_t ctx, slong prec)
{
  slong n = fmpq_mpoly_ctx_nvars(ctx);
  slong L = fmpq_mpoly_length(pol, ctx);
  slong* degrees = flint_malloc(n * sizeof(slong));
  slong j, k;
  acb_ptr* powers = flint_malloc(n * sizeof(acb_ptr));
  acb_t res, temp;
  fmpq_t coeff;
  slong exp;

  fmpq_mpoly_degrees_si(degrees, pol, ctx);
  for (k = 0; k < n; k++)
    {
      powers[k] = _acb_vec_init(degrees[k]+1);
      acb_one(&(powers[k][0]));
      for (j = 1; j <= degrees[k]; j++)
	{
	  acb_mul(&(powers[k][j]), &(powers[k][j-1]), &vals[k], prec);
	}
    }
  acb_init(res);
  acb_init(temp);
  fmpq_init(coeff);

  acb_zero(res);
  for (j = 0; j < L; j++)
    {
      fmpq_mpoly_get_term_coeff_fmpq(coeff, pol, j, ctx);
      acb_set_fmpq(temp, coeff, prec);
      for (k = 0; k < n; k++)
	{
	  exp = fmpq_mpoly_get_term_var_exp_si(pol, j, k, ctx);
	  acb_mul(temp, temp, &(powers[k][exp]), prec);
	}
      acb_add(res, res, temp, prec);
    }

  acb_set(ev, res);

  acb_clear(res);
  acb_clear(temp);
  fmpq_clear(coeff);
  for (k = 0; k < n; k++)
    {
      _acb_vec_clear(powers[k], degrees[k]+1);
    }
  flint_free(degrees);
  flint_free(powers);
}


#include "igusa.h"

void igusa_from_curve(acb_ptr I, const acb_poly_t crv, slong prec)
{
  acb_ptr ai;
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  acb_ptr S;
  
  ai = _acb_vec_init(7);
  fmpq_mpoly_ctx_init(ctx, 7, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(7);
  S = _acb_vec_init(4);

  hdme_data_vars_set(vars, "a", 0);
  hdme_data_vars_set(vars, "b", 1);
  hdme_data_vars_set(vars, "c", 2);
  hdme_data_vars_set(vars, "d", 3);
  hdme_data_vars_set(vars, "e", 4);
  hdme_data_vars_set(vars, "f", 5);
  hdme_data_vars_set(vars, "g", 6);

  curve_coeffs(ai, crv);

  hdme_data_read(pol, (const char**) vars, "igusa/I4", ctx);
  hdme_data_evaluate_acb(&S[0], pol, ai, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "igusa/I6prime", ctx);
  hdme_data_evaluate_acb(&S[1], pol, ai, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "igusa/I10", ctx);
  hdme_data_evaluate_acb(&S[2], pol, ai, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "igusa/I2", ctx);
  hdme_data_evaluate_acb(&S[3], pol, ai, ctx, prec);
  acb_mul(&S[3], &S[3], &S[2], prec);

  igusa_from_streng(I, S, prec);

  _acb_vec_clear(ai, 7);
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 7);
  _acb_vec_clear(S, 4);  
}

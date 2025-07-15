
#include "igusa.h"

void igusa_R2_from_IC(acb_t res, acb_srcptr I, slong prec)
{  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  acb_ptr vals;
  
  fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(4);
  vals = _acb_vec_init(4);

  _acb_vec_set(vals, I, 4);

  hdme_data_vars_set(vars, "a", 0);
  hdme_data_vars_set(vars, "b", 1);
  hdme_data_vars_set(vars, "c", 2);
  hdme_data_vars_set(vars, "d", 3);
  
  hdme_data_read(pol, (const char**) vars, "igusa/R2", ctx);
  hdme_data_evaluate_acb(res, pol, vals, ctx, prec);

  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 4);
  _acb_vec_clear(vals, 4);
}

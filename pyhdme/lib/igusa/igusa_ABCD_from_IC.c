
#include "igusa.h"

void igusa_ABCD_from_IC(acb_ptr ABCD, acb_srcptr I, slong prec)
{  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  acb_ptr vals;
  
  fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(4);
  vals = _acb_vec_init(4);

  hdme_data_vars_set(vars, "A", 0);
  hdme_data_vars_set(vars, "B", 1);
  hdme_data_vars_set(vars, "C", 2);
  hdme_data_vars_set(vars, "D", 3);

  _acb_vec_set(vals, I, 4);
    
  hdme_data_read(pol, (const char**) vars, "igusa/A", ctx);
  hdme_data_evaluate_acb(&vals[0], pol, vals, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "igusa/B", ctx);
  hdme_data_evaluate_acb(&vals[1], pol, vals, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "igusa/C", ctx);
  hdme_data_evaluate_acb(&vals[2], pol, vals, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "igusa/D", ctx);
  hdme_data_evaluate_acb(&vals[3], pol, vals, ctx, prec);
  
  _acb_vec_set(ABCD, vals, 4);

  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 4);
  _acb_vec_clear(vals, 4);
}

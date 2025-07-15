
#include "igusa.h"

void igusa_ABCD_from_IC_fmpz(fmpq* ABCD, fmpz* IC)
{  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  fmpq* vals;
  fmpz_t one;
  slong j;
  
  fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(4);
  vals = _fmpq_vec_init(4);
  fmpz_init(one);

  hdme_data_vars_set(vars, "A", 0);
  hdme_data_vars_set(vars, "B", 1);
  hdme_data_vars_set(vars, "C", 2);
  hdme_data_vars_set(vars, "D", 3);

  fmpz_one(one);
  for (j = 0; j < 4; j++) fmpq_set_fmpz_frac(&vals[j], &IC[j], one);
    
  hdme_data_read(pol, (const char**) vars, "igusa/A", ctx);
  hdme_data_evaluate_fmpq(&vals[0], pol, vals, ctx);
  hdme_data_read(pol, (const char**) vars, "igusa/B", ctx);
  hdme_data_evaluate_fmpq(&vals[1], pol, vals, ctx);
  hdme_data_read(pol, (const char**) vars, "igusa/C", ctx);
  hdme_data_evaluate_fmpq(&vals[2], pol, vals, ctx);
  hdme_data_read(pol, (const char**) vars, "igusa/D", ctx);
  hdme_data_evaluate_fmpq(&vals[3], pol, vals, ctx);
  
  for (j = 0; j < 4; j++) fmpq_set(&ABCD[j], &vals[j]);

  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 4);
  _fmpq_vec_clear(vals, 4);
  fmpz_clear(one);
}


#include "igusa.h"

static acb_struct* A11(acb_ptr A) {return &A[0];}
static acb_struct* A12(acb_ptr A) {return &A[1];}
static acb_struct* A22(acb_ptr A) {return &A[2];}
static acb_struct* A33(acb_ptr A) {return &A[3];}

void cardona_conic(acb_ptr Aij, acb_srcptr ABCD, slong prec)
{  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  acb_ptr res;
  
  fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(4);
  res = _acb_vec_init(4);
  
  hdme_data_vars_set(vars, "A", 0);
  hdme_data_vars_set(vars, "B", 1);
  hdme_data_vars_set(vars, "C", 2);
  hdme_data_vars_set(vars, "D", 3);
  
  hdme_data_read(pol, (const char**) vars, "cardona/A11", ctx);
  hdme_data_evaluate_acb(A11(res), pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/A12", ctx);
  hdme_data_evaluate_acb(A12(res), pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/A22", ctx);
  hdme_data_evaluate_acb(A22(res), pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/A33", ctx);
  hdme_data_evaluate_acb(A33(res), pol, ABCD, ctx, prec);

  _acb_vec_set(Aij, res, 4);

  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 4);
  _acb_vec_clear(res, 4);
}

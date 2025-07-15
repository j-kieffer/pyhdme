
#include "igusa.h"

static acb_struct* a111(acb_ptr a) {return &a[0];}
static acb_struct* a112(acb_ptr a) {return &a[1];}
static acb_struct* a122(acb_ptr a) {return &a[2];}
static acb_struct* a133(acb_ptr a) {return &a[3];}
static acb_struct* a222(acb_ptr a) {return &a[4];}
static acb_struct* a233(acb_ptr a) {return &a[5];}

void cardona_cubic(acb_ptr aijk, acb_srcptr ABCD, slong prec)
{  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  acb_ptr res;
  
  fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(4);
  res = _acb_vec_init(6);
  
  hdme_data_vars_set(vars, "A", 0);
  hdme_data_vars_set(vars, "B", 1);
  hdme_data_vars_set(vars, "C", 2);
  hdme_data_vars_set(vars, "D", 3);

  hdme_data_read(pol, (const char**) vars, "cardona/a111", ctx);
  hdme_data_evaluate_acb(a111(res), pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a112", ctx);
  hdme_data_evaluate_acb(a112(res), pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a122", ctx);
  hdme_data_evaluate_acb(a122(res), pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a133", ctx);
  hdme_data_evaluate_acb(a133(res), pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a222", ctx);
  hdme_data_evaluate_acb(a222(res), pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a233", ctx);
  hdme_data_evaluate_acb(a233(res), pol, ABCD, ctx, prec);

  _acb_vec_set(aijk, res, 6);
    
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 4);
  _acb_vec_clear(res, 6);
}


#include "hilbert.h"

void humbert_AA1BB1B2(acb_ptr AA1BB1B2, acb_srcptr rs, slong delta,
		      slong prec)
{
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;

  fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(2);
  
  humbert_vars_set(vars, delta);
  
  if (delta == 33)
    {
      humbert_get_mpoly(pol, (const char**) vars, "I2", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[0], pol, rs, ctx, prec);
      humbert_get_mpoly(pol, (const char**) vars, "I4", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[1], pol, rs, ctx, prec);
      humbert_get_mpoly(pol, (const char**) vars, "I6", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[2], pol, rs, ctx, prec);
      humbert_get_mpoly(pol, (const char**) vars, "I10", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[3], pol, rs, ctx, prec);
      humbert_get_mpoly(pol, (const char**) vars, "den", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[4], pol, rs, ctx, prec);
    }
  else
    {
      humbert_get_mpoly(pol, (const char**) vars, "A", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[0], pol, rs, ctx, prec);
      
      /* acb_printd(&rs[0], 30); flint_printf("\n"); */
      /* acb_printd(&rs[1], 30); flint_printf("\n"); */
      /* fmpq_mpoly_print_pretty(pol, (const char**) vars, ctx); flint_printf("\n"); */
      /* acb_printd(&AA1BB1B2[0], 30); flint_printf("\n"); */
      
      humbert_get_mpoly(pol, (const char**) vars, "A1", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[1], pol, rs, ctx, prec);      
      humbert_get_mpoly(pol, (const char**) vars, "B", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[2], pol, rs, ctx, prec);      
      humbert_get_mpoly(pol, (const char**) vars, "B1", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[3], pol, rs, ctx, prec);
      if (delta == 61)
	{   
	  humbert_get_mpoly(pol, (const char**) vars, "B1den", delta, ctx);
	  hdme_data_evaluate_acb(&AA1BB1B2[4], pol, rs, ctx, prec);
	  acb_div(&AA1BB1B2[3], &AA1BB1B2[3], &AA1BB1B2[4], prec);
	}    
      humbert_get_mpoly(pol, (const char**) vars, "B2", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[4], pol, rs, ctx, prec);
    }
  
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 2);
}

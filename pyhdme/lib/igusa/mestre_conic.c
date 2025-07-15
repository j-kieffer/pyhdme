
#include "igusa.h"

void mestre_conic(acb_ptr conic, acb_srcptr ABCD, const acb_t U, const acb_t I10, slong prec)
{
  acb_ptr A_coefs;
  acb_ptr C_coefs;
  acb_t temp;
  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;

  A_coefs = _acb_vec_init(6);
  C_coefs = _acb_vec_init(6);
  acb_init(temp);
  
  fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(4);

  hdme_data_vars_set(vars, "A", 0);
  hdme_data_vars_set(vars, "B", 1);
  hdme_data_vars_set(vars, "C", 2);
  hdme_data_vars_set(vars, "D", 3);
  
  hdme_data_read(pol, (const char**) vars, "mestre/conic/A11", ctx);
  hdme_data_evaluate_acb(&A_coefs[0], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/conic/A22", ctx);
  hdme_data_evaluate_acb(&A_coefs[1], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/conic/A33", ctx);
  hdme_data_evaluate_acb(&A_coefs[2], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/conic/A23", ctx);
  hdme_data_evaluate_acb(&A_coefs[3], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/conic/A31", ctx);
  hdme_data_evaluate_acb(&A_coefs[4], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/conic/A12", ctx);
  hdme_data_evaluate_acb(&A_coefs[5], pol, ABCD, ctx, prec);

  /* A_coefs is now computed; compute C_coefs by hand */

  _acb_vec_set(C_coefs, A_coefs, 6);
  
  acb_pow_ui(temp, U, 2, prec);
  acb_mul(&C_coefs[0], &C_coefs[0], temp, prec);
  acb_pow_ui(temp, I10, 8, prec);
  acb_mul(&C_coefs[0], &C_coefs[0], temp, prec);
  
  acb_pow_ui(temp, I10, 10, prec);
  acb_mul(&C_coefs[1], &C_coefs[1], temp, prec);
  
  acb_pow_ui(temp, U, 8, prec);
  acb_mul(&C_coefs[2], &C_coefs[2], temp, prec);
  
  acb_pow_ui(temp, U, 4, prec);
  acb_mul(&C_coefs[3], &C_coefs[3], temp, prec);
  acb_pow_ui(temp, I10, 5, prec);
  acb_mul(&C_coefs[3], &C_coefs[3], temp, prec);
  
  acb_pow_ui(temp, U, 5, prec);
  acb_mul(&C_coefs[4], &C_coefs[4], temp, prec);
  acb_pow_ui(temp, I10, 4, prec);
  acb_mul(&C_coefs[4], &C_coefs[4], temp, prec);
  
  acb_pow_ui(temp, U, 1, prec);
  acb_mul(&C_coefs[5], &C_coefs[5], temp, prec);
  acb_pow_ui(temp, I10, 9, prec);
  acb_mul(&C_coefs[5], &C_coefs[5], temp, prec);
  
  _acb_vec_set(conic, C_coefs, 6);
  
  _acb_vec_clear(A_coefs, 6);
  _acb_vec_clear(C_coefs, 6);
  acb_clear(temp);
  
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 4);
}


#include "igusa.h"

void mestre_cubic(acb_ptr cubic, acb_srcptr ABCD, const acb_t U, const acb_t I10, slong prec)
{
  acb_ptr A_coefs;
  acb_ptr C_coefs;
  acb_t temp;
  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;

  A_coefs = _acb_vec_init(10);
  C_coefs = _acb_vec_init(10);
  acb_init(temp);
  fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(4);

  hdme_data_vars_set(vars, "A", 0);
  hdme_data_vars_set(vars, "B", 1);
  hdme_data_vars_set(vars, "C", 2);
  hdme_data_vars_set(vars, "D", 3);
  
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c111", ctx);
  hdme_data_evaluate_acb(&A_coefs[0], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c112", ctx);
  hdme_data_evaluate_acb(&A_coefs[1], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c113", ctx);
  hdme_data_evaluate_acb(&A_coefs[2], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c122", ctx);
  hdme_data_evaluate_acb(&A_coefs[3], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c123", ctx);
  hdme_data_evaluate_acb(&A_coefs[4], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c133", ctx);
  hdme_data_evaluate_acb(&A_coefs[5], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c222", ctx);
  hdme_data_evaluate_acb(&A_coefs[6], pol, ABCD, ctx, prec);  
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c223", ctx);
  hdme_data_evaluate_acb(&A_coefs[7], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c233", ctx);
  hdme_data_evaluate_acb(&A_coefs[8], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "mestre/cubic/c333", ctx);
  hdme_data_evaluate_acb(&A_coefs[9], pol, ABCD, ctx, prec);

  /* A_coefs is now set; compute C_coefs */
  _acb_vec_set(C_coefs, A_coefs, 10);

  /* 111 */
  acb_pow_ui(temp, U, 3, prec);
  acb_mul(&C_coefs[0], &C_coefs[0], temp, prec);
  acb_pow_ui(temp, I10, 12, prec);
  acb_mul(&C_coefs[0], &C_coefs[0], temp, prec);

  /* 112 */
  acb_pow_ui(temp, U, 2, prec);
  acb_mul(&C_coefs[1], &C_coefs[1], temp, prec);
  acb_pow_ui(temp, I10, 13, prec);
  acb_mul(&C_coefs[1], &C_coefs[1], temp, prec);

  /* 113 */
  acb_pow_ui(temp, U, 6, prec);
  acb_mul(&C_coefs[2], &C_coefs[2], temp, prec);
  acb_pow_ui(temp, I10, 8, prec);
  acb_mul(&C_coefs[2], &C_coefs[2], temp, prec);

  /* 122 */
  acb_pow_ui(temp, U, 1, prec);
  acb_mul(&C_coefs[3], &C_coefs[3], temp, prec);
  acb_pow_ui(temp, I10, 14, prec);
  acb_mul(&C_coefs[3], &C_coefs[3], temp, prec);

  /* 123 */
  acb_pow_ui(temp, U, 5, prec);
  acb_mul(&C_coefs[4], &C_coefs[4], temp, prec);
  acb_pow_ui(temp, I10, 9, prec);
  acb_mul(&C_coefs[4], &C_coefs[4], temp, prec);

  /* 133 */
  acb_pow_ui(temp, U, 9, prec);
  acb_mul(&C_coefs[5], &C_coefs[5], temp, prec);
  acb_pow_ui(temp, I10, 4, prec);
  acb_mul(&C_coefs[5], &C_coefs[5], temp, prec);

  /* 222 */
  acb_pow_ui(temp, I10, 15, prec);
  acb_mul(&C_coefs[6], &C_coefs[6], temp, prec);

  /* 223 */
  acb_pow_ui(temp, U, 4, prec);
  acb_mul(&C_coefs[7], &C_coefs[7], temp, prec);
  acb_pow_ui(temp, I10, 10, prec);
  acb_mul(&C_coefs[7], &C_coefs[7], temp, prec);

  /* 233 */
  acb_pow_ui(temp, U, 8, prec);
  acb_mul(&C_coefs[8], &C_coefs[8], temp, prec);
  acb_pow_ui(temp, I10, 5, prec);
  acb_mul(&C_coefs[8], &C_coefs[8], temp, prec);

  /* 333 */
  acb_pow_ui(temp, U, 12, prec);
  acb_mul(&C_coefs[9], &C_coefs[9], temp, prec);

  _acb_vec_set(cubic, C_coefs, 10);
  _acb_vec_clear(A_coefs, 10);
  _acb_vec_clear(C_coefs, 10);
  acb_clear(temp);
  
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 4);
}

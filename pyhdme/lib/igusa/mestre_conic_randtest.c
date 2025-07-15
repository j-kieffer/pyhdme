
#include "igusa.h"

/* Obtain the conic as (nonzero poly of degree 1)^2 + nonzero constant */
/* Set x1 = 1, x2 = t, x3 = t^10 to use univariate polys only */

void mestre_conic_randtest(acb_ptr conic, flint_rand_t state, slong prec)
{
  acb_poly_t aux, res;
  acb_t coef;
  slong mag_bits = 2;
  slong kron = 10;

  acb_poly_init(aux);
  acb_poly_init(res);
  acb_init(coef);

  acb_poly_set_coeff_si(aux, 0, 1);
  acb_randtest_precise(coef, state, prec, mag_bits);
  while (acb_contains_zero(coef)) acb_randtest_precise(coef, state, prec, mag_bits);
  acb_poly_set_coeff_acb(aux, 1, coef);
  acb_randtest_precise(coef, state, prec, mag_bits);
  while (acb_contains_zero(coef)) acb_randtest_precise(coef, state, prec, mag_bits);
  acb_poly_set_coeff_acb(aux, kron, coef);
  acb_poly_mul(res, aux, aux, prec);
  
  acb_randtest_precise(coef, state, prec, mag_bits);
  while (acb_contains_zero(coef)) acb_randtest_precise(coef, state, prec, mag_bits);
  acb_poly_set_acb(aux, coef);
  acb_poly_add(res, res, aux, prec);

  /* Coeffs in conic are: A11, A22, A33, A23, A31, A12 */
  acb_poly_fit_length(res, 2*kron+1);
  acb_set(&conic[0], acb_poly_get_coeff_ptr(res, 0));
  acb_set(&conic[1], acb_poly_get_coeff_ptr(res, 2));
  acb_set(&conic[2], acb_poly_get_coeff_ptr(res, 2*kron));
  acb_div_si(&conic[3], acb_poly_get_coeff_ptr(res, kron+1), 2, prec);
  acb_div_si(&conic[4], acb_poly_get_coeff_ptr(res, kron), 2, prec);
  acb_div_si(&conic[5], acb_poly_get_coeff_ptr(res, 1), 2, prec);

  acb_poly_clear(aux);
  acb_poly_clear(res);
  acb_clear(coef);
}

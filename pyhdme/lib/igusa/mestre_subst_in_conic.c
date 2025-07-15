
#include "igusa.h"

void mestre_subst_in_conic(acb_poly_t subst, const acb_poly_t x1, const acb_poly_t x2,
			   const acb_poly_t x3, acb_srcptr conic, slong prec)
{
  
  acb_t coeff;
  acb_poly_t res;
  acb_poly_t aux;

  acb_poly_init(res);
  acb_poly_init(aux);
  acb_init(coeff);

  /* Conic is C11*x1^2+C22*x2^2+C33*x3^2+2*C12*x1*x2+2*C13*x1*x3+2*C23*x2*x3 */
  /* Coefficients are stored in the order: 11 22 33 23 31 12 */

  /* Diagonal coefficients */
  acb_set(coeff, &conic[0]);
  acb_poly_mul(aux, x1, x1, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);

  acb_poly_mul(aux, x2, x2, prec);
  acb_set(coeff, &conic[1]);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);

  acb_poly_mul(aux, x3, x3, prec);
  acb_set(coeff, &conic[2]);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);

  /* Off-diagonal coefficients */
  acb_poly_mul(aux, x2, x3, prec);
  acb_mul_si(coeff, &conic[3], 2, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  acb_poly_mul(aux, x3, x1, prec);
  acb_mul_si(coeff, &conic[4], 2, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  acb_poly_mul(aux, x1, x2, prec);
  acb_mul_si(coeff, &conic[5], 2, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);

  /* Set result */
  acb_poly_set(subst, res);
  
  acb_poly_clear(res);
  acb_poly_clear(aux);
  acb_clear(coeff);
}

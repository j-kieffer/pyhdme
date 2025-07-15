
#include "igusa.h"

/* Factor computations of x1sqr, x2sqr, x3sqr ? */
void mestre_subst_in_cubic(acb_poly_t subst, const acb_poly_t x1, const acb_poly_t x2,
			   const acb_poly_t x3, acb_srcptr cubic, slong prec)
{
  acb_t coeff;
  acb_poly_t res;
  acb_poly_t aux;

  acb_poly_init(res);
  acb_poly_init(aux);
  acb_init(coeff);

  /* Cubic is
     c111*x1^3+c222*x2^3+c333*x3^3+3*c112*x1^2*x2+3*c113*x1^2*x3+
     3*c122*x1*x2^2+3*c133*x1*x3^2+3*c233*x2*x3^2+3*c223*x2^2*x3+
     6*c123*x1*x2*x3 */
  /* Coefficients are stored in the order c111 c112 c113 c122 c123
     c133 c222 c223 c233 c333 */

  /* Diagonal coefficients */
  acb_set(coeff, &cubic[0]);
  acb_poly_pow_ui(aux, x1, 3, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  acb_set(coeff, &cubic[6]);
  acb_poly_pow_ui(aux, x2, 3, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  acb_set(coeff, &cubic[9]);
  acb_poly_pow_ui(aux, x3, 3, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);

  /* Crossed coefficients  */
  acb_mul_si(coeff, &cubic[1], 3, prec);
  acb_poly_mul(aux, x1, x1, prec);
  acb_poly_mul(aux, aux, x2, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  acb_mul_si(coeff, &cubic[2], 3, prec);
  acb_poly_mul(aux, x1, x1, prec);
  acb_poly_mul(aux, aux, x3, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  acb_mul_si(coeff, &cubic[3], 3, prec);
  acb_poly_mul(aux, x2, x2, prec);
  acb_poly_mul(aux, aux, x1, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  acb_mul_si(coeff, &cubic[5], 3, prec);
  acb_poly_mul(aux, x3, x3, prec);
  acb_poly_mul(aux, aux, x1, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  acb_mul_si(coeff, &cubic[7], 3, prec);
  acb_poly_mul(aux, x2, x2, prec);
  acb_poly_mul(aux, aux, x3, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  acb_mul_si(coeff, &cubic[8], 3, prec);
  acb_poly_mul(aux, x3, x3, prec);
  acb_poly_mul(aux, aux, x2, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);

  /* Mixed coefficient */
  acb_mul_si(coeff, &cubic[4], 6, prec);
  acb_poly_mul(aux, x1, x2, prec);
  acb_poly_mul(aux, aux, x3, prec);
  acb_poly_scalar_mul(aux, aux, coeff, prec);
  acb_poly_add(res, res, aux, prec);
  
  /* Set result */
  acb_poly_set(subst, res);
  
  acb_poly_clear(res);
  acb_poly_clear(aux);
  acb_clear(coeff);
}

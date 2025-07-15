
#include "igusa.h"

void mestre_eval_cubic(acb_t res, acb_srcptr pt, acb_srcptr cubic, slong prec)
{
  acb_t ev, term;

  /* Cubic is
     c111*x1^3+c222*x2^3+c333*x3^3+3*c112*x1^2*x2+3*c113*x1^2*x3+
     3*c122*x1*x2^2+3*c133*x1*x3^2+3*c233*x2*x3^2+3*c223*x2^2*x3+
     6*c123*x1*x2*x3 */
  /* Coefficients are stored in the order c111 c112 c113 c122 c123
     c133 c222 c223 c233 c333 */
  acb_init(ev);
  acb_init(term);

  acb_pow_si(term, &pt[0], 3, prec);
  acb_addmul(ev, term, &cubic[0], prec);

  acb_sqr(term, &pt[0], prec);
  acb_mul(term, term, &pt[1], prec);
  acb_mul_si(term, term, 3, prec);
  acb_addmul(ev, term, &cubic[1], prec);

  acb_sqr(term, &pt[0], prec);
  acb_mul(term, term, &pt[2], prec);
  acb_mul_si(term, term, 3, prec);
  acb_addmul(ev, term, &cubic[2], prec);

  acb_sqr(term, &pt[1], prec);
  acb_mul(term, term, &pt[0], prec);
  acb_mul_si(term, term, 3, prec);
  acb_addmul(ev, term, &cubic[3], prec);

  acb_mul(term, &pt[0], &pt[1], prec);
  acb_mul(term, term, &pt[2], prec);
  acb_mul_si(term, term, 6, prec);
  acb_addmul(ev, term, &cubic[4], prec);

  acb_sqr(term, &pt[2], prec);
  acb_mul(term, term, &pt[0], prec);
  acb_mul_si(term, term, 3, prec);
  acb_addmul(ev, term, &cubic[5], prec);

  acb_pow_si(term, &pt[1], 3, prec);
  acb_addmul(ev, term, &cubic[6], prec);

  acb_sqr(term, &pt[1], prec);
  acb_mul(term, term, &pt[2], prec);
  acb_mul_si(term, term, 3, prec);
  acb_addmul(ev, term, &cubic[7], prec);

  acb_sqr(term, &pt[2], prec);
  acb_mul(term, term, &pt[1], prec);
  acb_mul_si(term, term, 3, prec);
  acb_addmul(ev, term, &cubic[8], prec);
  
  acb_pow_si(term, &pt[2], 3, prec);
  acb_addmul(ev, term, &cubic[9], prec);
  
  acb_set(res, ev);
			   
  acb_clear(ev);
  acb_clear(term);
}

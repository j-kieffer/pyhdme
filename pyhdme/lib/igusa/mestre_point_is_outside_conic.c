
#include "igusa.h"

int mestre_point_is_outside_conic(acb_srcptr pt, acb_srcptr conic, slong prec)
{
  acb_t ev, term;
  int res = 0;

  /* Recall coefficients are stored in the following order:
     A11, A22, A33, A23, A31, A12
     and pt contains x1, x2, x3. */
  acb_init(ev);
  acb_init(term);

  acb_sqr(term, &pt[0], prec);
  acb_addmul(ev, term, &conic[0], prec);
  acb_sqr(term, &pt[1], prec);
  acb_addmul(ev, term, &conic[1], prec);
  acb_sqr(term, &pt[2], prec);
  acb_addmul(ev, term, &conic[2], prec);
  
  acb_mul(term, &pt[1], &pt[2], prec);
  acb_mul_si(term, term, 2, prec);
  acb_addmul(ev, term, &conic[3], prec);
  acb_mul(term, &pt[0], &pt[2], prec);
  acb_mul_si(term, term, 2, prec);
  acb_addmul(ev, term, &conic[4], prec);
  acb_mul(term, &pt[0], &pt[1], prec);
  acb_mul_si(term, term, 2, prec);
  acb_addmul(ev, term, &conic[5], prec);

  res = !acb_contains_zero(ev);
			   
  acb_clear(ev);
  acb_clear(term);
  return res;
}

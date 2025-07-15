
#include "igusa.h"

static void
igusa_I6prime_from_IC(acb_t I6prime, acb_srcptr I, slong prec)
{
  acb_t res, temp;

  acb_init(res);
  acb_init(temp);
  
  /* Get I6prime from I6 */
  acb_mul_si(res, &I[2], -3, prec);
  acb_mul(temp, &I[0], &I[1], prec);
  acb_add(res, res, temp, prec);
  acb_div_si(res, res, 2, prec);

  acb_set(I6prime, res);
  acb_clear(res);
  acb_clear(temp);
}


void igusa_from_IC(acb_ptr I, acb_srcptr IC, slong prec)
{
  acb_t I12, I6prime;
  acb_ptr S;

  acb_init(I12);
  acb_init(I6prime);
  S = _acb_vec_init(4);

  acb_mul(I12, &IC[0], &IC[3], prec);
  igusa_I6prime_from_IC(I6prime, IC, prec);

  acb_set(&S[0], &IC[1]);
  acb_set(&S[1], I6prime);
  acb_set(&S[2], &IC[3]);
  acb_set(&S[3], I12);

  igusa_from_streng(I, S, prec);

  acb_clear(I12);
  acb_clear(I6prime);
  _acb_vec_clear(S, 4);
}

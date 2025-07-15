
#include "igusa.h"

static void
igusa_I2_from_streng(acb_t I2, acb_srcptr S, slong prec)
{
  acb_div(I2, &S[3], &S[2], prec);
}

static void
igusa_I6_from_streng(acb_t I6, acb_srcptr S, slong prec)
{
  acb_t I2, res, temp;

  acb_init(I2);
  acb_init(res);
  acb_init(temp);

  acb_div(I2, &S[3], &S[2], prec);
  acb_mul_si(res, &S[1], 2, prec);
  acb_mul(temp, I2, &S[0], prec);
  acb_sub(res, res, temp, prec);
  acb_div_si(res, res, -3, prec);
  
  acb_set(I6, res);
  
  acb_clear(I2);
  acb_clear(res);
  acb_clear(temp);
}

void igusa_IC(acb_ptr IC, acb_srcptr I, slong prec)
{
  acb_t I2, I6;
  acb_ptr S;

  acb_init(I2);
  acb_init(I6);
  S = _acb_vec_init(4);

  igusa_streng(S, I, prec);

  igusa_I2_from_streng(I2, S, prec);
  igusa_I6_from_streng(I6, S, prec);

  acb_set(&IC[0], I2);
  acb_set(&IC[1], &S[0]);
  acb_set(&IC[2], I6);
  acb_set(&IC[3], &S[2]);

  acb_clear(I2);
  acb_clear(I6);
  _acb_vec_clear(S, 4);
}

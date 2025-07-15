
#include "hilbert.h"

void humbert_parametrize(acb_ptr I, acb_srcptr rs, slong delta,
			 slong prec)
{
  acb_ptr AA1BB1B2;
  acb_ptr res;
  acb_t temp;

  AA1BB1B2 = _acb_vec_init(5);
  res = _acb_vec_init(4);
  acb_init(temp);

  humbert_AA1BB1B2(AA1BB1B2, rs, delta, prec);
  if (delta == 33)
    {
      _acb_vec_set(res, AA1BB1B2, 4);
      acb_div(&res[0], &res[0], &AA1BB1B2[4], prec);
      acb_div(&res[2], &res[2], &AA1BB1B2[4], prec);
      acb_mul(temp, &res[0], &res[1], prec);
      acb_mul_si(&res[2], &res[2], -3, prec);
      acb_add(&res[2], &res[2], temp, prec);
      acb_div_si(&res[2], &res[2], 2, prec);
      
      acb_set(&I[0], &res[1]);
      acb_set(&I[1], &res[2]);
      acb_set(&I[2], &res[3]);
      acb_mul(&I[3], &res[3], &res[0], prec);  
      igusa_from_streng(I, I, prec);
    }
  else
    {
      humbert_igusa_from_AA1BB1B2(I, AA1BB1B2, prec);
    }

  _acb_vec_clear(AA1BB1B2, 5);
  _acb_vec_clear(res, 4);
  acb_clear(temp);
}

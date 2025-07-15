
#include "igusa.h"

void cov_rescale(acb_ptr I, acb_srcptr S, const acb_t scal,
		 slong nb, slong* weights, slong prec)
{
  acb_t temp;
  acb_ptr res;
  slong j;
  
  acb_init(temp);
  res = _acb_vec_init(nb);

  for (j = 0; j < nb; j++)
    {
      acb_pow_ui(temp, scal, weights[j], prec);
      acb_mul(&res[j], &S[j], temp, prec);
    }
  
  _acb_vec_set(I, res, nb);
  _acb_vec_clear(res, nb);
  acb_clear(temp);
}


#include "hilbert.h"

void hilbert_scalar_mul(acb_ptr z, const fmpz_poly_t x, acb_srcptr t,
			slong delta, slong prec)
{
  acb_ptr res;
  acb_t c;

  res = _acb_vec_init(2);
  acb_init(c);

  hilbert_sigma1(c, x, delta, prec);
  acb_mul(&res[0], &t[0], c, prec);
  hilbert_sigma2(c, x, delta, prec);
  acb_mul(&res[1], &t[1], c, prec);

  _acb_vec_set(z, res, 2);
  _acb_vec_clear(res, 2);
  acb_clear(c);
}

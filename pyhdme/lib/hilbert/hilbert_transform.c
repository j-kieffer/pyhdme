
#include "hilbert.h"

void hilbert_transform(acb_ptr z, const fmpz_poly_mat_t m, acb_srcptr t,
		       slong delta, slong prec)
{
  acb_t a, b, c, d, temp, sqrtd;
  acb_ptr res;
  
  acb_init(a);
  acb_init(b);
  acb_init(c);
  acb_init(d);
  acb_init(temp);
  acb_init(sqrtd);
  res = _acb_vec_init(2);

  hilbert_sigma1(a, fmpz_poly_mat_entry(m, 0, 0), delta, prec);
  hilbert_sigma1(b, fmpz_poly_mat_entry(m, 0, 1), delta, prec);
  hilbert_sigma1(c, fmpz_poly_mat_entry(m, 1, 0), delta, prec);
  hilbert_sigma1(d, fmpz_poly_mat_entry(m, 1, 1), delta, prec);
  arb_sqrt_ui(acb_realref(sqrtd), delta, prec);
  acb_div(b, b, sqrtd, prec);
  acb_mul(c, c, sqrtd, prec);
  
  acb_mul(&res[0], a, &t[0], prec);
  acb_add(&res[0], &res[0], b, prec);
  acb_mul(temp, c, &t[0], prec);
  acb_add(temp, temp, d, prec);
  acb_div(&res[0], &res[0], temp, prec);

  hilbert_sigma2(a, fmpz_poly_mat_entry(m, 0, 0), delta, prec);
  hilbert_sigma2(b, fmpz_poly_mat_entry(m, 0, 1), delta, prec);
  hilbert_sigma2(c, fmpz_poly_mat_entry(m, 1, 0), delta, prec);
  hilbert_sigma2(d, fmpz_poly_mat_entry(m, 1, 1), delta, prec);
  acb_neg(sqrtd, sqrtd);
  acb_div(b, b, sqrtd, prec);
  acb_mul(c, c, sqrtd, prec);
  
  acb_mul(&res[1], a, &t[1], prec);
  acb_add(&res[1], &res[1], b, prec);
  acb_mul(temp, c, &t[1], prec);
  acb_add(temp, temp, d, prec);
  acb_div(&res[1], &res[1], temp, prec);

  _acb_vec_set(z, res, 2);

  acb_clear(a);
  acb_clear(b);
  acb_clear(c);
  acb_clear(d);
  acb_clear(temp);
  acb_clear(sqrtd);
  _acb_vec_clear(res, 2);
}

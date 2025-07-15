
#include "hilbert.h"

void hilbert_map(acb_mat_t tau, acb_srcptr t, slong delta, slong prec)
{
  acb_mat_t R, res;
  acb_mat_init(R, 2, 2);
  acb_mat_init(res, 2, 2);

  hilbert_R(R, delta, prec);
  acb_set(acb_mat_entry(res, 0, 0), &t[0]);
  acb_set(acb_mat_entry(res, 1, 1), &t[1]);
  acb_mat_mul(res, res, R, prec);

  acb_mat_transpose(R, R);
  acb_mat_mul(res, R, res, prec);

  acb_mat_set(tau, res);

  acb_mat_clear(R);
  acb_mat_clear(res);
}


#include "igusa.h"

int tau_from_igusa(acb_mat_t tau, acb_srcptr I, slong prec)
{
  acb_ptr theta2;
  int res;

  theta2 = _acb_vec_init(16);
  res = tau_theta2_from_igusa(tau, theta2, I, prec);
  _acb_vec_clear(theta2, 16);

  return res;
}

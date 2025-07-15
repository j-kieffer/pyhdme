
#include "igusa.h"

int igusa_from_tau(acb_ptr I, const acb_mat_t tau, slong prec)
{
  acb_ptr theta2;
  int res;

  theta2 = _acb_vec_init(16);
  
  res = theta2_unif(theta2, tau, prec);
  if (res) res = theta2_renormalize(theta2, theta2, prec);
  if (res) igusa_from_theta2(I, theta2, prec);
  
  _acb_vec_clear(theta2, 16);
  return res;
}

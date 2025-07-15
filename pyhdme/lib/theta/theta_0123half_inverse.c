
#include "theta.h"

int
theta_0123half_inverse(acb_mat_t tau, acb_srcptr th_half, slong prec)
{
  acb_ptr th;
  int res;
  
  th = _acb_vec_init(16);
  
  theta_duplication(th, th_half, prec);
  res = theta2_inverse(tau, th, prec);

  _acb_vec_clear(th, 16);
  return res;
}

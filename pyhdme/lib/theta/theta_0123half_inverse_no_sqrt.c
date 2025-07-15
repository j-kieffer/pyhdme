
#include "theta.h"

/* Use theta2_inverse_no_sqrt instead. */

int
theta_0123half_inverse_no_sqrt(acb_mat_t tau, acb_srcptr th_half, slong prec)
{
  acb_ptr th;
  int res;
  
  th = _acb_vec_init(16);
  
  theta_duplication(th, th_half, prec);
  res = theta2_inverse_no_sqrt(tau, th, prec);

  _acb_vec_clear(th, 16);
  return res;
}

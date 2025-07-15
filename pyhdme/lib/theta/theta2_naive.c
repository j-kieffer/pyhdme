
#include "theta.h"

int
theta2_naive(acb_ptr th, const acb_mat_t tau, slong prec)
{
  acb_mat_t tau_half;
  acb_ptr th_half;

  int res;

  acb_mat_init(tau_half, 2, 2);
  th_half = _acb_vec_init(4);
  
  acb_mat_scalar_mul_2exp_si(tau_half, tau, -1);
  res = theta_0123_naive(th_half, tau_half, prec);
  theta_duplication(th, th_half, prec);

  acb_mat_clear(tau_half);
  _acb_vec_clear(th_half, 4);

  return res;
}

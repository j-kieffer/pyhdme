
#include "theta.h"

int theta2_newton(acb_ptr th2, const acb_mat_t tau, slong prec)
{
  acb_ptr current_th;
  acb_ptr next_th;
  acb_mat_t tau_half;
  acb_t th0;
  mag_t err;
  
  slong g = 2;
  slong delta = THETA_NEWTON_LOSS;
  slong current_prec, next_prec;
  slong j;
  int res = 1;
  int v = get_theta_verbose();

  acb_mat_init(tau_half, g, g);
  current_th = _acb_vec_init(4);
  next_th = _acb_vec_init(4);
  acb_init(th0);
  mag_init(err);

  current_prec = theta2_newton_start_prec(prec);
  acb_mat_scalar_mul_2exp_si(tau_half, tau, -1);
  theta_0123_naive(current_th, tau_half, current_prec);
  /* Normalize so that first value is 1 */
  acb_set(th0, &current_th[0]);
  _acb_vec_scalar_div(current_th, current_th, 4, th0, prec);

  if (v > 0) {flint_printf("(theta2_newton) %wd: ", prec); fflush(stdout);}
  
  while ((current_prec < prec) && res)
    {
      if (v > 0) {flint_printf("%wd, ", current_prec); fflush(stdout);}
      next_prec = 2*current_prec;
      for (j = 0; j < 4; j++) acb_get_mid(&next_th[j], &current_th[j]);
      res = theta2_newton_step(next_th, tau, next_th, next_prec);
      /* next_th has actually next_prec - delta correct bits */
      current_prec = next_prec - delta;
      _acb_vec_set(current_th, next_th, 4);
    }
  if (v > 0) {flint_printf("%wd.\n", current_prec); fflush(stdout);}
  
  /* Add some heuristic error */
  mag_set_ui_2exp_si(err, 1, -prec+delta);
  for (j = 0; j < 4; j++) acb_add_error_mag(&current_th[j], err);
  /* current_th contains th_0123half at the desired precision, use duplication */
  theta_duplication(th2, current_th, prec);

  acb_mat_clear(tau_half);
  _acb_vec_clear(current_th, 4);
  _acb_vec_clear(next_th, 4);
  acb_clear(th0);
  mag_clear(err);
  return res;
}

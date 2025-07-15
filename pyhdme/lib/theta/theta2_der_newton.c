
#include "theta.h"

int theta2_der_newton(acb_ptr th2, acb_mat_t dth2, const acb_mat_t tau, slong prec)
{
  acb_ptr current_th;
  acb_ptr next_th;
  acb_mat_t tau_half;
  acb_mat_t dth_approx;
  acb_mat_t dth_4x3;
  acb_t th0;
  mag_t err;
  
  slong g = 2;
  slong delta = THETA_NEWTON_LOSS;
  slong current_prec, next_prec;
  slong j, k;
  int res = 1;
  int v = get_theta_verbose();

  acb_mat_init(tau_half, g, g);
  current_th = _acb_vec_init(4);
  next_th = _acb_vec_init(4);
  acb_mat_init(dth_approx, 3, 3);
  acb_mat_init(dth_4x3, 4, 3);
  acb_init(th0);
  mag_init(err);

  current_prec = theta2_newton_start_prec(prec);
  acb_mat_scalar_mul_2exp_si(tau_half, tau, -1);
  theta_0123_naive(current_th, tau_half, current_prec);
  /* Normalize so that first value is 1 */
  acb_set(th0, &current_th[0]);
  _acb_vec_scalar_div(current_th, current_th, 4, th0, prec);

  if (v > 0) {flint_printf("(theta2_der_newton) %wd: ", prec); fflush(stdout);}
  
  while ((current_prec < prec) && res)
    {
      if (v > 0) {flint_printf("%wd, ", current_prec); fflush(stdout);}
      next_prec = 2*current_prec;
      for (j = 0; j < 4; j++) acb_get_mid(&next_th[j], &current_th[j]);
      res = theta2_der_newton_step(next_th, dth_approx, tau, next_th, next_prec);
      /* next_th has actually next_prec - delta correct bits */
      current_prec = next_prec - delta;
      _acb_vec_set(current_th, next_th, 4);
    }
  if (v > 0) {flint_printf("%wd.\n", current_prec); fflush(stdout);}
  
  /* Add some heuristic error */
  mag_set_ui_2exp_si(err, 1, -prec+delta);
  for (j = 0; j < 4; j++) acb_add_error_mag(&current_th[j], err);
  /* current_th contains th_0123half at the desired precision */
  /* dth_approx contains th_0123half_diff to the desired precision,
     with first coordinate fixed to be 1; use duplication */
  
  if (acb_mat_is_zero(dth_approx))
    {
      theta2_der_naive(th2, dth2, tau, prec);
    }
  else
    {
      mag_set_ui_2exp_si(err, 1, -prec/2+delta);
      for (j = 0; j < 3; j++)
	{
	  for (k = 0; k < 3; k++)
	    {
	      acb_mul_si(acb_mat_entry(dth_4x3, j+1, k),
			 acb_mat_entry(dth_approx, j, k), 2, prec);
	      acb_add_error_mag(acb_mat_entry(dth_4x3, j+1, k), err);
	    }
	  acb_mul(acb_mat_entry(dth_4x3, j+1, 2),
		  acb_mat_entry(dth_4x3, j+1, 2),
		  acb_mat_entry(tau, 0, 1), prec);
	  acb_mul_si(acb_mat_entry(dth_4x3, j+1, 2),
		     acb_mat_entry(dth_4x3, j+1, 2), 2, prec);
	}
      theta_der_duplication(th2, dth2, current_th, dth_4x3, prec);
    }

  acb_mat_clear(tau_half);
  _acb_vec_clear(current_th, 4);
  _acb_vec_clear(next_th, 4);
  acb_mat_clear(dth_approx);
  acb_mat_clear(dth_4x3);
  acb_clear(th0);
  mag_clear(err);
  return res;
}

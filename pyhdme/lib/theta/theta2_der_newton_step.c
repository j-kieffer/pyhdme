
#include "theta.h"

int theta2_der_newton_step(acb_ptr th_half, acb_mat_t dth_approx,
			   const acb_mat_t tau, acb_srcptr th_half_approx,
			   slong prec)
{
  acb_mat_t dtau;
  acb_mat_t tau_approx;
  acb_mat_t dth_half;
  acb_ptr th_half_corr;
  int res;
  slong i;

  acb_mat_init(dtau, 3, 3);
  acb_mat_init(dth_half, 3, 3);
  acb_mat_init(tau_approx, 2, 2);
  th_half_corr = _acb_vec_init(3);

  /* Compute tau_approx corresponding to th_half_approx */
  res = theta_0123half_inverse(tau_approx, th_half_approx, prec);
  if (!res)
    {
      flint_printf("(theta2_der_newton_step) Warning: unable to invert the following theta values\n");
      for (i = 0; i < 3; i++)
	{
	  flint_printf("    "); acb_printd(&th_half_approx[i], 30); flint_printf("\n");
	}
    }

  
  /* Compute dtau at these approx values using finite differences */
  if (res)
    {
      res = theta_0123half_inverse_diff(dtau, tau_approx, th_half_approx, prec);
    }
  /* Invert matrix to get dth_half at tau_approx */
  if (res)
    {
      res = acb_mat_inv(dth_half, dtau, prec);
    }
  
  /* Compute corrections to apply to th_half */
  acb_sub(acb_mat_entry(tau_approx, 0, 0), acb_mat_entry(tau_approx, 0, 0),
	  acb_mat_entry(tau, 0, 0), prec);
  acb_sub(acb_mat_entry(tau_approx, 1, 1), acb_mat_entry(tau_approx, 1, 1),
	  acb_mat_entry(tau, 1, 1), prec);
  acb_sqr(acb_mat_entry(tau_approx, 1, 0), acb_mat_entry(tau_approx, 1, 0), prec);
  acb_submul(acb_mat_entry(tau_approx, 1, 0), acb_mat_entry(tau, 1, 0),
	     acb_mat_entry(tau, 1, 0), prec);
  
  for (i = 0; i < 3; i++)
    {
      acb_mul(&th_half_corr[i],
	      acb_mat_entry(dth_half, i, 0),
	      acb_mat_entry(tau_approx, 0, 0), prec);
      acb_addmul(&th_half_corr[i],
		 acb_mat_entry(dth_half, i, 1),
		 acb_mat_entry(tau_approx, 1, 1), prec);
      acb_addmul(&th_half_corr[i],
		 acb_mat_entry(dth_half, i, 2),
		 acb_mat_entry(tau_approx, 1, 0), prec);
      acb_neg(&th_half_corr[i], &th_half_corr[i]);
    }
  /* Apply corrections */
  _acb_vec_set(th_half, th_half_approx, 4);
  for (i = 0; i < 3; i++)
    {
      acb_add(&th_half[i+1], &th_half[i+1], &th_half_corr[i], prec);
    }
  acb_mat_set(dth_approx, dth_half);
  
  acb_mat_clear(dtau);
  acb_mat_clear(dth_half);
  acb_mat_clear(tau_approx);
  _acb_vec_clear(th_half_corr, 3);
  return res;
}

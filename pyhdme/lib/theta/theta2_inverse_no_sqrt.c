
#include "theta.h"

/* Same as theta2_inverse, BUT we leave tau_3^2 as is. */

int
theta2_inverse_no_sqrt(acb_mat_t tau, acb_srcptr th, slong prec)
{
  acb_ptr means;
  acb_ptr a;
  acb_ptr tau_entries;
  
  int i, j;
  int res = 1;

  int theta_labels[4][4] = {
    {0, 1, 2, 3},
    {4, 0, 6, 2},
    {8, 9, 0, 1},
    {0, 8, 4, 12}
  };

  means = _acb_vec_init(4);
  a = _acb_vec_init(4);
  tau_entries = _acb_vec_init(3);

  /* Compute 1/theta_0^2(gamma_i tau) */
  for (i = 0; i < 4; i++)
    {
      /* Set means[i] to a Borchardt mean */
      for (j = 0; j < 4; j++)
	{
	  acb_set(&a[j], &th[theta_char_set_label_g2(theta_labels[i][j])]);
	}
      res = res && borchardt_mean(&means[i], a, prec);
    }

  if (res)
    {
      /* tau1 = i m[0]/m[1], tau2 = i m[0]/m[2], tau3^2 - tau1*tau2 = m[0]/m[3] */
      acb_div(&tau_entries[0], &means[0], &means[1], prec);
      acb_mul_onei(&tau_entries[0], &tau_entries[0]);
      acb_div(&tau_entries[1], &means[0], &means[2], prec);
      acb_mul_onei(&tau_entries[1], &tau_entries[1]);
      acb_div(&tau_entries[2], &means[0], &means[3], prec);
      acb_addmul(&tau_entries[2], &tau_entries[0], &tau_entries[1], prec);
    }
  
  /* Set tau */
  acb_set(acb_mat_entry(tau, 0, 0), &tau_entries[0]);
  acb_set(acb_mat_entry(tau, 1, 1), &tau_entries[1]);
  acb_set(acb_mat_entry(tau, 0, 1), &tau_entries[2]);
  acb_set(acb_mat_entry(tau, 1, 0), &tau_entries[2]);
  
  _acb_vec_clear(means, 4);
  _acb_vec_clear(a, 4);
  _acb_vec_clear(tau_entries, 3);
  return res;
}

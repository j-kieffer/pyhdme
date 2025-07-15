
#include "theta.h"

int
theta2_inverse(acb_mat_t tau, acb_srcptr th, slong prec)
{
  acb_ptr means;
  acb_ptr a;
  acb_ptr tau_entries;
  arb_t im;
  
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
  arb_init(im);

  /* Compute 1/theta_0^2(gamma_i tau) */
  for (i = 0; i < 4; i++)
    {
      /* Set means[i] to a Borchardt mean */
      for (j = 0; j < 4; j++)
	{
	  acb_set(&a[j], &th[theta_char_set_label_g2(theta_labels[i][j])]);
	}
      res = res && borchardt_mean(&means[i], a, prec);
      /* Debug */
      /* flint_printf("(theta2_inverse) Values:\n");
      for (j = 0; j < 4; j++)
	{
	  flint_printf("    "); acb_printd(&a[j], 30); flint_printf("\n");
	}
      flint_printf("(theta2_inverse) Borchardt mean:");
      acb_printd(&means[i], 30);
      flint_printf("\n"); */
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

      /* Debug */
      /* flint_printf("(theta2_inverse) Value before sqrt:\n    ");
      acb_printd(&tau_entries[2], 30);
      flint_printf("\n"); */
      
      borchardt_sqrt(&tau_entries[2], &tau_entries[2], prec);
      /* Adjust the sign of tau3 */
      arb_set(im, acb_imagref(&tau_entries[2]));
      if (arf_cmp_si(arb_midref(im), 0) < 0)
	{
	  acb_neg(&tau_entries[2], &tau_entries[2]);
	}
    }
  
  /* Set tau */
  acb_set(acb_mat_entry(tau, 0, 0), &tau_entries[0]);
  acb_set(acb_mat_entry(tau, 1, 1), &tau_entries[1]);
  acb_set(acb_mat_entry(tau, 0, 1), &tau_entries[2]);
  acb_set(acb_mat_entry(tau, 1, 0), &tau_entries[2]);

  /* Debug */
  /* flint_printf("(theta2_inverse) Result:\n");
     acb_mat_printd(tau, 10); */
  
  _acb_vec_clear(means, 4);
  _acb_vec_clear(a, 4);
  _acb_vec_clear(tau_entries, 3);
  arb_clear(im);
  return res;
}

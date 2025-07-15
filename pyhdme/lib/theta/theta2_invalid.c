
#include "theta.h"

int theta2_invalid(acb_srcptr th2, slong prec)
{
  acb_ptr means;
  acb_ptr a;
  acb_ptr tau_entries;
  acb_mat_t tau;
  arb_t im;
  
  int res = 0;
  int b_success = 1;
  slong i, j;
  
  int theta_labels[4][4] = {
    {0, 1, 2, 3},
    {4, 0, 6, 2},
    {8, 9, 0, 1},
    {0, 8, 4, 12}
  };
  
  means = _acb_vec_init(4);
  a = _acb_vec_init(4);
  tau_entries = _acb_vec_init(3);
  acb_mat_init(tau, 2, 2);
  arb_init(im);

  /* Either one of the four Borchardt means is definitely invalid; or
     the Borchardt means succeed and the resulting matrix does not
     belong to the fundamental domain. */
  for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 4; j++)
	{
	  acb_set(&a[j], &th2[theta_char_set_label_g2(theta_labels[i][j])]);
	}
      res = res || borchardt_mean_invalid(a, prec);
      if (!res && b_success)
	{
	  b_success = borchardt_mean(&means[i], a, prec);
	  /* if (!b_success)
	    {
	      flint_printf("(theta2_invalid) Undecided Borchardt mean:\n");
	      for (j = 0; j < 4; j++)
		{
		  acb_printd(&a[j], 30); flint_printf("\n");
		}
		} */
	}
    }
  
  if (!res && b_success)
    {
      acb_div(&tau_entries[0], &means[0], &means[1], prec);
      acb_mul_onei(&tau_entries[0], &tau_entries[0]);
      acb_div(&tau_entries[1], &means[0], &means[2], prec);
      acb_mul_onei(&tau_entries[1], &tau_entries[1]);
      acb_div(&tau_entries[2], &means[0], &means[3], prec);
      acb_addmul(&tau_entries[2], &tau_entries[0], &tau_entries[1], prec);
      acb_sqrt(&tau_entries[2], &tau_entries[2], prec);

      arb_set(im, acb_imagref(&tau_entries[2]));
      if (arf_cmp_si(arb_midref(im), 0) < 0)
	{
	  acb_neg(&tau_entries[2], &tau_entries[2]);
	}
      
      acb_set(acb_mat_entry(tau, 0, 0), &tau_entries[0]);
      acb_set(acb_mat_entry(tau, 1, 1), &tau_entries[1]);
      acb_set(acb_mat_entry(tau, 0, 1), &tau_entries[2]);
      acb_set(acb_mat_entry(tau, 1, 0), &tau_entries[2]);

      /* acb_mat_printd(tau, 10); */
      
      res = siegel_not_in_fundamental_domain(tau, prec);
    }
  
  _acb_vec_clear(means, 4);
  _acb_vec_clear(a, 4);
  _acb_vec_clear(tau_entries, 3);
  acb_mat_clear(tau);
  arb_clear(im);
  return res;
}


#include "theta.h"


void theta_der_duplication(acb_ptr th2_2tau, acb_mat_t th2_der_2tau,
			   acb_srcptr th_tau, const acb_mat_t th_der_tau,
			   slong prec)
{
  ulong a, b, b1, b2, ch, ch1, ch2;
  acb_t thb1, thb2;
  acb_ptr dthb1;
  acb_ptr res;
  acb_mat_t dres;
  slong g = 2;
  slong k;
  
  acb_init(thb1);
  acb_init(thb2);
  dthb1 = _acb_vec_init(3);
  res = _acb_vec_init(n_pow(2, 2*g));
  acb_mat_init(dres, n_pow(2, 2*g), 3);
  
  /* Compute each multiplication only once */
  for (b = 0; b < n_pow(2, g); b++)
    {
      for (b1 = 0; b1 < n_pow(2, g); b1++)
	{
	  b2 = b ^ b1; /* Addition mod 2 */
	  ch1 = theta_char_set_ab(0, b1, g);
	  acb_set(thb1, &th_tau[ch1]);
	  ch2 = theta_char_set_ab(0, b2, g);
	  acb_set(thb2, &th_tau[ch2]);
	  
	  acb_mul(thb1, thb1, thb2, prec);
	  acb_mul_2exp_si(thb1, thb1, -g);

	  /* Compute derivative of thb1 wrt entries tau */
	  for (k = 0; k < 3; k++)
	    {
	      acb_mul(&dthb1[k], &th_tau[ch1],
		      acb_mat_entry(th_der_tau, ch2, k), prec);
	      acb_addmul(&dthb1[k],
			 acb_mat_entry(th_der_tau, ch1, k), &th_tau[ch2], prec);
	      acb_mul_2exp_si(&dthb1[k], &dthb1[k], -g);
	    }
	  
	  for (a = 0; a < n_pow(2, g); a++)
	    {
	      ch = theta_char_set_ab(a, b, g);
	      if (!theta_char_is_even(ch, g)) continue;
	      
	      if (theta_char_dot_product(a, b1, g) % 2 == 0)
		{
		  acb_add(&res[ch], &res[ch], thb1, prec);
		  for (k = 0; k < 3; k++)
		    {
		      acb_add(acb_mat_entry(dres, ch, k),
			      acb_mat_entry(dres, ch, k), &dthb1[k], prec);
		    }
		}
	      else
		{
		  acb_sub(&res[ch], &res[ch], thb1, prec);
		  for (k = 0; k < 3; k++)
		    {
		      acb_sub(acb_mat_entry(dres, ch, k),
			      acb_mat_entry(dres, ch, k), &dthb1[k], prec);
		    }
		}
	    }
	}
    }
  /* We want derivatives wrt entries of 2tau */
  acb_mat_scalar_div_si(dres, dres, 2, prec);
  acb_mat_set(th2_der_2tau, dres);
  _acb_vec_set(th2_2tau, res, n_pow(2, 2*g));
  
  acb_clear(thb1);
  acb_clear(thb2);
  _acb_vec_clear(dthb1, 3);
  _acb_vec_clear(res, n_pow(2, 2*g));
  acb_mat_clear(dres);
}


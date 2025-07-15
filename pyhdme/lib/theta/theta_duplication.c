
#include "theta.h"

void
theta_duplication(acb_ptr th2_2tau, acb_srcptr th_tau, slong prec)
{
  ulong a, b, b1, b2, ch;
  acb_t thb1, thb2;
  acb_ptr res;
  slong g = 2;
  
  acb_init(thb1);
  acb_init(thb2);
  res = _acb_vec_init(n_pow(2, 2*g));
  
  /* Compute each multiplication only once */
  for (b = 0; b < n_pow(2, g); b++)
    {
      for (b1 = 0; b1 < n_pow(2, g); b1++)
	{
	  b2 = b ^ b1; /* Addition mod 2 */
	  ch = theta_char_set_ab(0, b1, g);
	  acb_set(thb1, &th_tau[ch]);
	  ch = theta_char_set_ab(0, b2, g);
	  acb_set(thb2, &th_tau[ch]);
	  
	  acb_mul(thb1, thb1, thb2, prec);
	  acb_mul_2exp_si(thb1, thb1, -g);
	  
	  for (a = 0; a < n_pow(2, g); a++)
	    {
	      ch = theta_char_set_ab(a, b, g);
	      if (!theta_char_is_even(ch, g)) continue;
	      
	      if (theta_char_dot_product(a, b1, g) % 2 == 0)
		{
		  acb_add(&res[ch], &res[ch], thb1, prec);
		}
	      else
		{
		  acb_sub(&res[ch], &res[ch], thb1, prec);
		}
	    }
	}
    }
  _acb_vec_set(th2_2tau, res, n_pow(2, 2*g));
  
  acb_clear(thb1);
  acb_clear(thb2);
  _acb_vec_clear(res, n_pow(2, 2*g));
}

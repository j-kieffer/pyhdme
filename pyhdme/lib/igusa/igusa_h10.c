
#include "igusa.h"

void igusa_h10(acb_t h10, acb_srcptr theta2, slong prec)
{ 
  slong g = 2;
  ulong ch;
  ulong n = n_pow(2, 2*g);

  acb_t res;
  
  acb_init(res);
  acb_one(res);

  for (ch = 0; ch < n; ch++)
    {
      if (theta_char_is_even(ch, g))
	{
	  acb_mul(res, res, &theta2[ch], prec);
	}
    }
  
  acb_set(h10, res);
  acb_clear(res);
}


#include "igusa.h"

void igusa_h4(acb_t h4, acb_srcptr theta2, slong prec)
{
  slong g = 2;
  ulong ch;
  acb_t res, aux;

  acb_init(res);
  acb_init(aux);

  for (ch = 0; ch < n_pow(2, 2*g); ch++)
    {
      if (theta_char_is_even(ch, g))
	{
	  acb_pow_ui(aux, &theta2[ch], 4, prec);
	  acb_add(res, res, aux, prec);
	}
    }

  acb_set(h4, res);
  acb_clear(res);
  acb_clear(aux);
}

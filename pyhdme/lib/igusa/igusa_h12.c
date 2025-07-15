
#include "igusa.h"

void igusa_h12(acb_t h12, acb_srcptr theta2, slong prec)
{
  slong g = 2;
  ulong ch1, ch2, ch3, ch4, ch;
  ulong n = n_pow(2, 2*g);

  acb_t res, aux;
  
  acb_init(res);
  acb_init(aux);

  for (ch1 = 0; ch1 < n; ch1++) {
    for (ch2 = ch1 + 1; ch2 < n; ch2++) {
      for (ch3 = ch2 + 1; ch3 < n; ch3++) {
	for (ch4 = ch3 + 1; ch4 < n; ch4++) {
	  if (theta_char_is_goepel(ch1, ch2, ch3, ch4, g)) {
	    acb_one(aux);
	    for (ch = 0; ch < n; ch++) {
	      if (theta_char_is_even(ch, g)
		  && (ch != ch1)
		  && (ch != ch2)
		  && (ch != ch3)
		  && (ch != ch4)) {
		acb_mul(aux, aux, &theta2[ch], prec);
	      }
	    }
	    acb_sqr(aux, aux, prec);
	    acb_add(res, res, aux, prec);
	  }
	}
      }
    }
  }
  
  acb_set(h12, res);
  acb_clear(res);
  acb_clear(aux);
}

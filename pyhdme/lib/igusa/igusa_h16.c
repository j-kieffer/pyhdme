
#include "igusa.h"

void igusa_h16(acb_t h16, acb_srcptr theta2, slong prec)
{
  slong g = 2;
  ulong ch1, ch2, ch3, ch4, ch;
  ulong n = n_pow(2, 2*g);

  acb_t res, aux1, aux2;
  
  acb_init(res);
  acb_init(aux1);
  acb_init(aux2);

  /* Compute aux1 as in igusa_h12; aux2 is the sum of theta^4 on the
     goepel quadruple */
  
  for (ch1 = 0; ch1 < n; ch1++) {
    for (ch2 = ch1 + 1; ch2 < n; ch2++) {
      for (ch3 = ch2 + 1; ch3 < n; ch3++) {
	for (ch4 = ch3 + 1; ch4 < n; ch4++) {
	  if (theta_char_is_goepel(ch1, ch2, ch3, ch4, g)) {
	    /* Start with aux2, using aux1 as temp */
	    acb_zero(aux2);
	    acb_pow_ui(aux1, &theta2[ch1], 4, prec);
	    acb_add(aux2, aux2, aux1, prec);
	    acb_pow_ui(aux1, &theta2[ch2], 4, prec);
	    acb_add(aux2, aux2, aux1, prec);
	    acb_pow_ui(aux1, &theta2[ch3], 4, prec);
	    acb_add(aux2, aux2, aux1, prec);
	    acb_pow_ui(aux1, &theta2[ch4], 4, prec);
	    acb_add(aux2, aux2, aux1, prec);
	    /* aux2 is done; compute aux1 */
	    acb_one(aux1);
	    for (ch = 0; ch < n; ch++) {
	      if (theta_char_is_even(ch, g)
		  && (ch != ch1)
		  && (ch != ch2)
		  && (ch != ch3)
		  && (ch != ch4)) {
		acb_mul(aux1, aux1, &theta2[ch], prec);
	      }
	    }
	    acb_sqr(aux1, aux1, prec);
	    /* Add in res */
	    acb_addmul(res, aux1, aux2, prec);
	  }
	}
      }
    }
  }
  
  acb_set(h16, res);
  acb_clear(res);
  acb_clear(aux1);
  acb_clear(aux2);
}


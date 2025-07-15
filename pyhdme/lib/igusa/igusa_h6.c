
#include "igusa.h"

static slong
igusa_h6_sign(ulong b, ulong c, ulong d)
{
  slong sgn;
  int b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3; /* d4 not used */

  b4 = b % 2; b = b >> 1;
  b3 = b % 2; b = b >> 1;
  b2 = b % 2; b = b >> 1;
  b1 = b % 2; b = b >> 1;
  
  c4 = c % 2; c = c >> 1;
  c3 = c % 2; c = c >> 1;
  c2 = c % 2; c = c >> 1;
  c1 = c % 2; c = c >> 1;
  
  d = d >> 1;
  d3 = d % 2; d = d >> 1;
  d2 = d % 2; d = d >> 1;
  d1 = d % 2; d = d >> 1;

  sgn = b1 + b2 + c1 + c2 + d1 + d2 + b1*c1 + b2*c2 + b4*c2 + b1*c3 - b2*c4 +
      b1*d1 - b3*d1 + c1*d1 + b2*d2 + c2*d2 + c4*d2 + c1*d3 - b2*b3*c1 -
      b2*b4*c2 - b1*b2*c3 - b2*b3*d1 - b3*c1*d1 - b1*c3*d1 - b2*c3*d1
      - b2*b4*d2 - b4*c2*d2 - b1*b2*d3 - b1*c1*d3 - b2*c1*d3;

  if ((sgn % 2) == 1) sgn = -1; else sgn = 1;
  
  return sgn;
}


void igusa_h6(acb_t h6, acb_srcptr theta2, slong prec)
{
  slong g = 2;
  ulong ch1, ch2, ch3;
  ulong n = n_pow(2, 2*g);
  slong sgn;

  acb_t res, aux;

  acb_init(res);
  acb_init(aux);

  for (ch1 = 0; ch1 < n; ch1++) {
    for (ch2 = ch1 + 1; ch2 < n; ch2++) {
      for (ch3 = ch2 + 1; ch3 < n; ch3++) {
	if (theta_char_is_syzygous(ch1, ch2, ch3, g))
	  {
	    sgn = igusa_h6_sign(ch1, ch2, ch3);
	    acb_mul(aux, &theta2[ch1], &theta2[ch2], prec);
	    acb_mul(aux, aux, &theta2[ch3], prec);
	    acb_sqr(aux, aux, prec);
	    acb_mul_si(aux, aux, sgn, prec);
	    acb_add(res, res, aux, prec);
	  }
      }
    }	
  }

  acb_set(h6, res);
  acb_clear(res);
  acb_clear(aux);
}

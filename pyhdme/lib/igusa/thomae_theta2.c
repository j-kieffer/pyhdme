
#include "igusa.h"

/* Factor square roots? */

/* We reduce to choosing 4 signs instead of 10, then apply Thomae's formulae */

void thomae_theta2(acb_ptr th2, acb_srcptr th4, acb_srcptr ros, slong signs, slong prec)
{
  slong k;
  slong ch;
  acb_t l, m, n;
  acb_t num, den, aux;
  acb_ptr res;

  acb_init(num);
  acb_init(den);
  acb_init(aux);
  acb_init(l);
  acb_init(m);
  acb_init(n);
  res = _acb_vec_init(16);

  acb_set(l, &ros[0]);
  acb_set(m, &ros[1]);
  acb_set(n, &ros[2]);

  acb_one(&res[0]);

  /* "Fundamental" ones (why?): 1, 2, 4, 8 */
  for (k = 1; k < 9; k *= 2) {
    acb_sqrt(&res[k], &th4[k], prec); /* th4[0] is one. */
    if (signs % 2 == 1) acb_neg(&res[k], &res[k]);
    signs = signs/2;
  }

  /* Set other entries using Thomae's formulae */
  ch = theta_char_set_label_g2(6);
  acb_set(num, &res[theta_char_set_label_g2(2)]);
  acb_set(den, n);
  acb_mul(den, den, &res[theta_char_set_label_g2(4)], prec);
  acb_div(&res[ch], num, den, prec);

  ch = theta_char_set_label_g2(12);
  acb_set(num, &res[theta_char_set_label_g2(8)]);
  acb_set(den, l);
  acb_mul(den, den, &res[theta_char_set_label_g2(4)], prec);
  acb_div(&res[ch], num, den, prec);

  ch = theta_char_set_label_g2(3);
  acb_mul(num, &res[theta_char_set_label_g2(4)],
	  &res[theta_char_set_label_g2(6)], prec);
  acb_sub_si(aux, n, 1, prec);
  acb_mul(num, num, aux, prec);
  acb_set(den, &res[theta_char_set_label_g2(1)]);
  acb_div(&res[ch], num, den, prec);

  ch = theta_char_set_label_g2(9);
  acb_mul(num, &res[theta_char_set_label_g2(4)],
	  &res[theta_char_set_label_g2(12)], prec);
  acb_sub_si(aux, l, 1, prec);
  acb_mul(num, num, aux, prec);
  acb_set(den, &res[theta_char_set_label_g2(1)]);
  acb_div(&res[ch], num, den, prec);

  ch = theta_char_set_label_g2(15);
  acb_mul(num, &res[theta_char_set_label_g2(1)],
	  &res[theta_char_set_label_g2(12)], prec);
  acb_sub(aux, n, m, prec);
  acb_mul(num, num, aux, prec);
  acb_set(den, &res[theta_char_set_label_g2(2)]);
  acb_sub_si(aux, n, 1, prec);
  acb_mul(den, den, aux, prec);
  acb_div(&res[ch], num, den, prec);

  _acb_vec_set(th2, res, 16);

  acb_clear(num);
  acb_clear(den);
  acb_clear(aux);
  acb_clear(l);
  acb_clear(m);
  acb_clear(n);
  _acb_vec_clear(res, 16);
}


#include "theta.h"
#include "igusa.h"

void thomae_theta4(acb_ptr th4, acb_srcptr ros, slong prec)
{
  acb_t l, m, n; /* Rosenhain invariants */
  acb_t num, den, aux;
  slong ch;
  
  acb_init(l);
  acb_init(m);
  acb_init(n);
  acb_init(num);
  acb_init(den);
  acb_init(aux);

  acb_set(l, &ros[0]);
  acb_set(m, &ros[1]);
  acb_set(n, &ros[2]);

  /* Thomae's formulas */
  ch = theta_char_set_label_g2(0);
  acb_one(&th4[ch]);

  ch = theta_char_set_label_g2(1);
  acb_set(num, m);
  acb_sub_si(aux, l, 1, prec);
  acb_mul(num, num, aux, prec);
  acb_sub_si(aux, n, 1, prec);
  acb_mul(num, num, aux, prec);
  acb_mul(den, l, n, prec);
  acb_sub_si(aux, m, 1, prec);
  acb_mul(den, den, aux, prec);
  acb_div(&th4[ch], num, den, prec);

  ch = theta_char_set_label_g2(2);
  acb_set(num, m);
  acb_sub_si(aux, l, 1, prec);
  acb_mul(num, num, aux, prec);
  acb_sub(aux, n, m, prec);
  acb_mul(num, num, aux, prec);
  acb_set(den, l);
  acb_sub_si(aux, m, 1, prec);
  acb_mul(den, den, aux, prec);
  acb_sub(aux, n, l, prec);
  acb_mul(den, den, aux, prec);
  acb_div(&th4[ch], num, den, prec);

  ch = theta_char_set_label_g2(4);
  acb_set(num, m);
  acb_mul(den, l, n, prec);
  acb_div(&th4[ch], num, den, prec);

  ch = theta_char_set_label_g2(8);
  acb_set(num, m);
  acb_sub(aux, l, m, prec);
  acb_mul(num, num, aux, prec);
  acb_sub_si(aux, n, 1, prec);
  acb_mul(num, num, aux, prec);
  acb_set(den, n);
  acb_sub_si(aux, m, 1, prec);
  acb_mul(den, den, aux, prec);
  acb_sub(aux, l, n, prec);
  acb_mul(den, den, aux, prec);
  acb_div(&th4[ch], num, den, prec);

  ch = theta_char_set_label_g2(6);
  acb_set(num, &th4[theta_char_set_label_g2(2)]);
  acb_sqr(den, n, prec);
  acb_mul(den, den, &th4[theta_char_set_label_g2(4)], prec);
  acb_div(&th4[ch], num, den, prec);

  ch = theta_char_set_label_g2(12);
  acb_set(num, &th4[theta_char_set_label_g2(8)]);
  acb_sqr(den, l, prec);
  acb_mul(den, den, &th4[theta_char_set_label_g2(4)], prec);
  acb_div(&th4[ch], num, den, prec);

  ch = theta_char_set_label_g2(3);
  acb_mul(num, &th4[theta_char_set_label_g2(4)],
	  &th4[theta_char_set_label_g2(6)], prec);
  acb_sub_si(aux, n, 1, prec);
  acb_sqr(aux, aux, prec);
  acb_mul(num, num, aux, prec);
  acb_set(den, &th4[theta_char_set_label_g2(1)]);
  acb_div(&th4[ch], num, den, prec);

  ch = theta_char_set_label_g2(9);
  acb_mul(num, &th4[theta_char_set_label_g2(4)],
	  &th4[theta_char_set_label_g2(12)], prec);
  acb_sub_si(aux, l, 1, prec);
  acb_sqr(aux, aux, prec);
  acb_mul(num, num, aux, prec);
  acb_set(den, &th4[theta_char_set_label_g2(1)]);
  acb_div(&th4[ch], num, den, prec);
  
  ch = theta_char_set_label_g2(15);
  acb_mul(num, &th4[theta_char_set_label_g2(1)],
	  &th4[theta_char_set_label_g2(12)], prec);
  acb_sub(aux, n, m, prec);
  acb_sqr(aux, aux, prec);
  acb_mul(num, num, aux, prec);
  acb_set(den, &th4[theta_char_set_label_g2(2)]);
  acb_sub_si(aux, n, 1, prec);
  acb_sqr(aux, aux, prec);
  acb_mul(den, den, aux, prec);
  acb_div(&th4[ch], num, den, prec);
  
  acb_clear(l);
  acb_clear(m);
  acb_clear(n);
  acb_clear(num);
  acb_clear(den);
  acb_clear(aux);
}

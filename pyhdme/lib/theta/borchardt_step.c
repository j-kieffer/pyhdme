
#include "theta.h"

/* Can we use only three square roots? */
int borchardt_step(acb_ptr a, acb_srcptr b, slong prec)
{
  int i;
  int res = 1;
  acb_ptr ra;

  ra = _acb_vec_init(4);

  /* Set a to b and compute square roots */
  for (i = 0; i < 4; i++)
    {
      acb_set(&a[i], &b[i]);
      res = res && acb_sqrt_goodpos(&ra[i], &a[i], prec);
    }

  if (res) /* Good square roots */
    {
      acb_add(&a[0], &a[0], &a[1], prec);
      acb_add(&a[0], &a[0], &a[2], prec);
      acb_add(&a[0], &a[0], &a[3], prec);
      acb_mul_2exp_si(&a[0], &a[0], -2);   /* a0 <- (a0+a1+a2+a3)/4 */

      acb_mul(&a[1], &ra[0], &ra[1], prec);
      acb_addmul(&a[1], &ra[2], &ra[3], prec);
      acb_mul_2exp_si(&a[1], &a[1], -1); /* a1 <- (\sqrt{a0}\sqrt{a1}+\sqrt{a2}\sqrt{a3})/2 */
      
      acb_mul(&a[2], &ra[0], &ra[2], prec);
      acb_addmul(&a[2], &ra[1], &ra[3], prec);
      acb_mul_2exp_si(&a[2], &a[2], -1);
      
      acb_mul(&a[3], &ra[0], &ra[3], prec);
      acb_addmul(&a[3], &ra[2], &ra[1], prec);
      acb_mul_2exp_si(&a[3], &a[3], -1);
    }

  _acb_vec_clear(ra, 4);
  return res;
}



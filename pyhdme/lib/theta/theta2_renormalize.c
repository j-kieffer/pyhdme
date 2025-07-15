
#include "theta.h"

int theta2_renormalize(acb_ptr th2, acb_srcptr th2_proj, slong prec)
{
  acb_ptr a;
  acb_t b;
  slong k;
  int res;

  a = _acb_vec_init(16);
  acb_init(b);

  _acb_vec_set(a, th2_proj, 16);
  for (k = 1; k < 16; k++)
    {
      acb_div(&a[k], &a[k], &a[0], prec);
    }
  acb_one(&a[0]);
  res = borchardt_mean(b, a, prec);
  if (res)
    {
      for (k = 0; k < 16; k++)
	{
	  acb_div(&a[k], &a[k], b, prec);
	}
    }
  _acb_vec_set(th2, a, 16);
  
  _acb_vec_clear(a, 16);
  acb_clear(b);
  return res;
}

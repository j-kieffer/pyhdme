
#include "hecke.h"

void hecke_slash(acb_ptr im, const acb_mat_t star, acb_srcptr val,
		 slong k, slong j, slong prec)
{
  acb_t det;

  acb_init(det);
  
  if (j != 0)
    {
      flint_printf("(hecke_slash) Error: not implemented for vector-valued forms (j = %wd)\n",
		   j);
      fflush(stdout);
      flint_abort();
    }
  else
    {
      acb_mat_det(det, star, prec);
      hecke_slash_scalar(im, det, val, k, prec);
    }
  
  acb_clear(det);
}

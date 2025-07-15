
#include "polynomials.h"

int acb_vec_contains_int(acb_srcptr v, slong nb)
{
  int res = 1;
  slong k;
  
  for (k = 0; k < nb; k++)
    {
      if (!acb_contains_int(&v[k]))
	{
	  res = 0;
	  break;
	}
    }

  return res;
}


#include "igusa.h"

static int cov_not_one(acb_srcptr I, slong nb)
{
  slong j;
  int res = 0;
  
  for (j = 0; j < nb; j++)
    {
      if (!arb_contains_si(acb_imagref(&I[j]), 0)
	  || !arb_contains_si(acb_realref(&I[j]), 1))
	{
	  res = 1;
	  break;
	}
    }
  
  return res;
}

int cov_no_rescale_to_one(acb_srcptr I, slong nb, slong* weights, slong prec)
{
  acb_t rt, zeta;
  acb_ptr test;
  slong j, k;
  slong i0 = -1;
  int res = 1;

  acb_init(rt);
  acb_init(zeta);
  test = _acb_vec_init(nb);

  /* Compute nonzero index */
  for (j = 0; j < nb; j++)
    {
      if (!acb_contains_zero(&I[j]))
	{
	  i0 = j;
	  break;
	}
    }
  if (i0 == -1) res = 0; /* Don't know. */

  if (res)
    {
      res = 1;
      /* Compute possible rescalings */
      borchardt_root_ui(rt, &I[i0], weights[i0], prec);
      acb_inv(rt, rt, prec);
      acb_unit_root(zeta, weights[i0], prec);
      for (k = 0; k < weights[i0]; k++)
	{
	  cov_rescale(test, I, rt, nb, weights, prec);
	  if (!cov_not_one(test, nb)) /* Found a potentially correct rescaling factor */
	    {
	      res = 0;
	      break;
	    }
	  acb_mul(rt, rt, zeta, prec);
	}
    }

  acb_clear(rt);
  acb_clear(zeta);
  _acb_vec_clear(test, nb);
  return res;
}

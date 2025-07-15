
#include "igusa.h"

int thomae_keep_candidate(const acb_mat_t tau, acb_srcptr I, slong prec)
{
  acb_ptr test;
  slong weights[4] = IGUSA_HALFWEIGHTS;
  int aux, res;

  test = _acb_vec_init(4);

  if (siegel_not_in_fundamental_domain(tau, prec))
    {
      res = 0;
    }
  else
    {
      aux = igusa_from_tau(test, tau, prec);
      res = 1;
      if (aux && cov_distinct(test, I, 4, weights, prec)) res = 0;
      /* We can say nothing if computing test fails */
    }
  
  _acb_vec_clear(test, 4);
  return res;
}

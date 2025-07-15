
#include "igusa.h"

void igusa_from_theta2(acb_ptr I, acb_srcptr theta2, slong prec)
{
  /* Support weird aliasing... */
  acb_ptr h;
  
  h = _acb_vec_init(4);
  
  igusa_h4(&h[0], theta2, prec);
  igusa_h6(&h[1], theta2, prec);
  igusa_h10(&h[2], theta2, prec);
  igusa_h12(&h[3], theta2, prec);

  acb_div_si(igusa_psi4(I), &h[0], 4, prec);
  acb_div_si(igusa_psi6(I), &h[1], 4, prec);
  acb_div_si(igusa_chi10(I), &h[2], -n_pow(2,12), prec);
  acb_div_si(igusa_chi12(I), &h[3], n_pow(2,15), prec);
  
  _acb_vec_clear(h, 4);
}


#include "igusa.h"

int tau_theta2_from_igusa(acb_mat_t tau, acb_ptr theta2, acb_srcptr I, slong prec)
{ 
  acb_ptr IC;
  acb_poly_t crv;
  int res;

  IC = _acb_vec_init(4);
  acb_poly_init(crv);

  igusa_IC(IC, I, prec);
  res = mestre(crv, IC, prec);
  if (res) res = tau_theta2_from_curve(tau, theta2, crv, prec);
  
  _acb_vec_clear(IC, 4);
  acb_poly_clear(crv);
  return res;
}

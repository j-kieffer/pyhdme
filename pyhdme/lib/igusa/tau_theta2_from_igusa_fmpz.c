
#include "igusa.h"

int tau_theta2_from_igusa_fmpz(acb_mat_t tau, acb_ptr th2, fmpz* I, slong prec)
{
  acb_poly_t crv;
  fmpz* IC;
  int res;

  IC = _fmpz_vec_init(4);
  acb_poly_init(crv);

  if (igusa_is_g2_curve_fmpz(I))
    {
      igusa_IC_fmpz(IC, I);
      res = mestre_fmpz(crv, IC, prec);
      if (res) res = tau_theta2_from_curve(tau, th2, crv, prec);
    }
  else
    {
      res = tau_theta2_from_igusa_ec(tau, th2, I, prec);
    }
  
  acb_poly_clear(crv);
  _fmpz_vec_clear(IC, 4);
  return res;  
}

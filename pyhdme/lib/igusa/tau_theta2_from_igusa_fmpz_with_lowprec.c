
#include "igusa.h"

/* See also tau_theta2_from_igusa_fmpz */

int tau_theta2_from_igusa_fmpz_with_lowprec(acb_mat_t tau, acb_ptr th2,
    fmpz* I, acb_srcptr th2_lp, slong prec)
{
  acb_poly_t crv;
  fmpz* IC;
  int res;

  IC = _fmpz_vec_init(4);
  acb_poly_init(crv);

  if (igusa_is_g2_curve_fmpz(I)) {
    igusa_IC_fmpz(IC, I);
    res = mestre_fmpz(crv, IC, prec);
    if (res)
      res = tau_theta2_from_curve_with_lowprec(tau, th2, crv, th2_lp, prec);
  } else {
    res = tau_theta2_from_igusa_ec(tau, th2, I, prec);
  }

  acb_poly_clear(crv);
  _fmpz_vec_clear(IC, 4);
  return res;
}

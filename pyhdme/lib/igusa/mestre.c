
#include "igusa.h"


int mestre(acb_poly_t crv, acb_srcptr IC, slong prec)
{
  acb_ptr conic;
  acb_ptr cubic;
  acb_ptr ABCD;
  acb_ptr pt;
  acb_poly_t x1, x2, x3;
  acb_t U;
  acb_t I10;
  int res = 1;

  
  if (!igusa_has_generic_automorphisms(IC, prec))
    {
      flint_printf("(mestre) Warning: cannot guarantee that curve does not have extra automorphisms\n");
      return 0;
    }

  conic = _acb_vec_init(6);
  cubic = _acb_vec_init(10);
  ABCD = _acb_vec_init(4);
  pt = _acb_vec_init(3);
  acb_poly_init(x1);
  acb_poly_init(x2);
  acb_poly_init(x3);
  acb_init(U);
  acb_init(I10);

  igusa_ABCD_from_IC(ABCD, IC, prec);
  acb_set(I10, &IC[3]);
  if (!acb_contains_zero(&ABCD[0])) acb_pow_si(U, &ABCD[0], 6, prec);
  else if (!acb_contains_zero(&ABCD[1])) acb_pow_si(U, &ABCD[1], 3, prec);
  else if (!acb_contains_zero(&ABCD[2])) acb_pow_si(U, &ABCD[2], 2, prec);
  else res = 0; /* This should never happen due to has_generic_automorphisms */

  if (res)
    {
      mestre_conic(conic, ABCD, U, I10, prec);
      mestre_cubic(cubic, ABCD, U, I10, prec);
      res = mestre_point_on_conic(pt, conic, prec);
    }

  if (res) /* Point on conic has been found */
    {
      mestre_parametrize_conic(x1, x2, x3, pt, conic, prec);
      mestre_subst_in_cubic(crv, x1, x2, x3, cubic, prec); 
    }

  if (res
      && !acb_contains_zero(acb_poly_get_coeff_ptr(crv, 0))
      && !acb_contains_zero(acb_poly_get_coeff_ptr(crv, 6))
      ) /* Balance coefficients: use U and x1 as temp */
    {
      acb_poly_get_coeff_acb(U, crv, 0);
      acb_poly_scalar_div(crv, crv, U, prec);
      borchardt_root_ui(U, acb_poly_get_coeff_ptr(crv, 6), 6, prec);
      acb_inv(U, U, prec);
      acb_poly_zero(x1);
      acb_poly_set_coeff_acb(x1, 1, U);
      acb_poly_compose(crv, crv, x1, prec);
    }

  _acb_vec_clear(conic, 6);
  _acb_vec_clear(cubic, 10);
  _acb_vec_clear(ABCD, 4);
  _acb_vec_clear(pt, 3);
  acb_poly_clear(x1);
  acb_poly_clear(x2);
  acb_poly_clear(x3);
  acb_clear(U);
  acb_clear(I10);
  return res;
}

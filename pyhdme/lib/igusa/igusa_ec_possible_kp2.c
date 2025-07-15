
#include "igusa.h"

int igusa_ec_possible_kp2(acb_ptr kp2, const acb_t j, slong prec)
{
  acb_poly_t pol;
  acb_poly_t term;
  acb_t c;
  slong nb, res;

  acb_poly_init(pol);
  acb_poly_init(term);
  acb_init(c);

  /* Set term to x^2(1-x)^2 */
  acb_poly_zero(term);
  acb_poly_set_coeff_si(term, 1, -1);
  acb_poly_set_coeff_si(term, 0, 1);
  acb_poly_mul(term, term, term, prec);
  acb_poly_shift_left(term, term, 2);

  acb_poly_scalar_mul(pol, term, j, prec);

  /* Set term to 256(1-x+x^2)^3 */
  acb_poly_zero(term);
  acb_poly_set_coeff_si(term, 2, 1);
  acb_poly_set_coeff_si(term, 1, -1);
  acb_poly_set_coeff_si(term, 0, 1);
  acb_poly_pow_ui(term, term, 3, prec);
  acb_set_si(c, 256);
  acb_poly_scalar_mul(term, term, c, prec);

  acb_poly_sub(pol, pol, term, prec);
  
  /* Has simple roots if j!=0, 1728 */
  /* See also thomae_roots */
  nb = acb_poly_find_roots(kp2, pol, NULL, prec, prec);
  res = (nb == 6);
  
  if (!res)
    {
      flint_printf("(igusa_possible_kp2) Warning: unable to isolate roots\n");
      acb_poly_printd(pol, 10); flint_printf("\n");
      acb_printd(j, 10); flint_printf("\n");
    }

  acb_poly_clear(pol);
  acb_poly_clear(term);
  acb_clear(c);
  return res;
}

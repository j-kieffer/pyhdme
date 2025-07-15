
#include "igusa.h"

/* We plug in 5 different values for a2, a3, b2, b3 - one of these is going to work. */
/* We don't want (0,0,0), so choose b2 and b3 nonzero. */
int mestre_point_on_conic(acb_ptr pt, acb_srcptr conic, slong prec)
{
  acb_poly_t subst_poly;
  acb_poly_t x1;
  acb_poly_t x2;
  acb_poly_t x3;
  acb_t c0, c1, c2;
  acb_t discr;
  acb_t lead;
  acb_ptr xx;
  
  slong a2, b2, a3, b3;
  int res = 0;

  acb_poly_init(subst_poly);
  acb_poly_init(x1);
  acb_poly_init(x2);
  acb_poly_init(x3);
  acb_init(c0);
  acb_init(c1);
  acb_init(c2);
  acb_init(discr);
  acb_init(lead);
  xx = _acb_vec_init(3);

  acb_poly_set_coeff_si(x1, 1, 1);
  
  for (a2 = 0; a2 < 5; a2++) {
    for (a3 = 0; a3 < 5; a3++) {
      for (b2 = 1; b2 < 6; b2++) {
	for (b3 = 1; b3 < 6; b3++) {
	  acb_poly_zero(x2);
	  acb_poly_zero(x3);
	  acb_poly_set_coeff_si(x2, 0, b2);
	  acb_poly_set_coeff_si(x2, 1, a2);
	  acb_poly_set_coeff_si(x3, 0, b3);
	  acb_poly_set_coeff_si(x3, 1, a3);
	  mestre_subst_in_conic(subst_poly, x1, x2, x3, conic, prec);
	  
	  /* Are discriminant and leading term nonzero ? */
	  acb_poly_get_coeff_acb(c0, subst_poly, 0);
	  acb_poly_get_coeff_acb(c1, subst_poly, 1);
	  acb_poly_get_coeff_acb(c2, subst_poly, 2);
	  acb_sqr(discr, c1, prec);
	  acb_mul(lead, c0, c2, prec); /* Use lead as temp */
	  acb_submul_si(discr, lead, 4, prec);
	  acb_poly_get_coeff_acb(lead, subst_poly, 2);
	  
	  if (!acb_contains_zero(discr)
	      && !acb_contains_zero(lead))
	    {
	      res = 1;
	      /* Compute x1, x2, x3 */
	      borchardt_sqrt(discr, discr, prec);
	      acb_sub(&xx[0], discr, c1, prec);
	      acb_div(&xx[0], &xx[0], lead, prec);
	      acb_div_si(&xx[0], &xx[0], 2, prec);
	      acb_mul_si(&xx[1], &xx[0], a2, prec);
	      acb_add_si(&xx[1], &xx[1], b2, prec);
	      acb_mul_si(&xx[2], &xx[0], a3, prec);
	      acb_add_si(&xx[2], &xx[2], b3, prec);
	      goto exit; /* Break all loops */
	    }
	}
      }
    }
  }
  goto exit;
  
 exit:
  {
    _acb_vec_set(pt, xx, 3);
    
    acb_poly_clear(subst_poly);
    acb_poly_clear(x1);
    acb_poly_clear(x2);
    acb_poly_clear(x3);
    acb_clear(c0);
    acb_clear(c1);
    acb_clear(c2);
    acb_clear(discr);
    acb_clear(lead);
    _acb_vec_clear(xx, 3);
    return res;
  }
}

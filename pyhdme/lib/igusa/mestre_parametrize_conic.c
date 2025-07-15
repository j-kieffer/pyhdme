
#include "igusa.h"

void mestre_parametrize_conic(acb_poly_t x1, acb_poly_t x2, acb_poly_t x3,
			     acb_srcptr pt, acb_srcptr conic, slong prec)
{
  slong i, j, k;
  acb_poly_t ui, uj, uk, aux;
  acb_poly_t subst;
  acb_poly_t a, b;
  slong kron = 10;

  acb_poly_init(ui);
  acb_poly_init(uj);
  acb_poly_init(uk);
  acb_poly_init(subst);
  acb_poly_init(aux);
  acb_poly_init(a);
  acb_poly_init(b);
  
  /* Find a nonzero coordinate in pt, called i */
  if (!acb_contains_zero(&pt[0])) {i = 0; j = 1; k = 2;}
  else if (!acb_contains_zero(&pt[1])) {i = 1; j = 0; k = 2;}
  else {i = 2; j = 0; k = 1;}

  /* Substitute in conic equation */
  /* Kind of Kronecker substitution: put t = m^kron */
  acb_poly_set_acb(ui, &pt[i]);
  acb_poly_set_coeff_acb(uj, 0, &pt[j]);
  acb_poly_set_coeff_si(uj, kron + 1, 1);
  acb_poly_set_coeff_acb(uk, 0, &pt[k]);
  acb_poly_set_coeff_si(uk, kron, 1);

  if (i == 0) mestre_subst_in_conic(subst, ui, uj, uk, conic, prec);
  else if (i == 1) mestre_subst_in_conic(subst, uj, ui, uk, conic, prec);
  else mestre_subst_in_conic(subst, uj, uk, ui, conic, prec);

  /* a = coeff of degree 2 in t (poly of degree 2 in m); 
     b = coeff of degree 1 in t (poly of degree 1 in m) */
  acb_poly_fit_length(subst, 2*kron + 3);
  acb_poly_set_coeff_acb(a, 0, acb_poly_get_coeff_ptr(subst, 2*kron));
  acb_poly_set_coeff_acb(a, 1, acb_poly_get_coeff_ptr(subst, 2*kron+1));
  acb_poly_set_coeff_acb(a, 2, acb_poly_get_coeff_ptr(subst, 2*kron+2));
  acb_poly_set_coeff_acb(b, 0, acb_poly_get_coeff_ptr(subst, kron));
  acb_poly_set_coeff_acb(b, 1, acb_poly_get_coeff_ptr(subst, kron+1));

  /* Get parametrization: xi, xj, xk */
  acb_poly_scalar_mul(ui, a, &pt[i], prec);
  acb_poly_scalar_mul(uj, a, &pt[j], prec);
  acb_poly_shift_left(aux, b, 1);
  acb_poly_sub(uj, uj, aux, prec);
  acb_poly_scalar_mul(uk, a, &pt[k], prec);
  acb_poly_sub(uk, uk, b, prec);

  /* Set result */
  if (i == 0)
    {
      acb_poly_set(x1, ui);
      acb_poly_set(x2, uj);
      acb_poly_set(x3, uk);
    }
  else if (i == 1)
    {
      acb_poly_set(x2, ui);
      acb_poly_set(x1, uj);
      acb_poly_set(x3, uk);
    }
  else
    {
      acb_poly_set(x3, ui);
      acb_poly_set(x1, uj);
      acb_poly_set(x2, uk);
    }
  
  acb_poly_clear(ui);
  acb_poly_clear(uj);
  acb_poly_clear(uk);
  acb_poly_clear(subst);
  acb_poly_clear(aux);
  acb_poly_clear(a);
  acb_poly_clear(b);
}

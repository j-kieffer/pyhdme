
#include "igusa.h"

static acb_struct* A11(acb_ptr A) {return &A[0];}
static acb_struct* A12(acb_ptr A) {return &A[1];}
static acb_struct* A22(acb_ptr A) {return &A[2];}
static acb_struct* A33(acb_ptr A) {return &A[3];}

static acb_struct* a111(acb_ptr a) {return &a[0];}
static acb_struct* a112(acb_ptr a) {return &a[1];}
static acb_struct* a122(acb_ptr a) {return &a[2];}
static acb_struct* a133(acb_ptr a) {return &a[3];}
static acb_struct* a222(acb_ptr a) {return &a[4];}
static acb_struct* a233(acb_ptr a) {return &a[5];}

void cardona(acb_poly_t crv, acb_srcptr IC, slong prec)
{
  acb_ptr ABCD;
  acb_ptr Aij; /* 11, 12, 22, 33 */
  acb_ptr aijk; /* 111, 112, 122, 133, 222, 233 */
  acb_poly_t P1, P2, P3;
  acb_t c;
  acb_poly_t term;
  
  ABCD = _acb_vec_init(4);
  Aij = _acb_vec_init(4);
  aijk = _acb_vec_init(6);
  acb_poly_init(P1);
  acb_poly_init(P2);
  acb_poly_init(P3);
  acb_init(c);
  acb_poly_init(term);

  /* Set coefficients */
  igusa_ABCD_from_IC(ABCD, IC, prec);
  cardona_conic(Aij, ABCD, prec);
  cardona_cubic(aijk, ABCD, prec);
  
  /* Set P1, P2, P3 */
  acb_mul_si(c, A12(Aij), -2, prec);
  acb_poly_set_coeff_acb(P1, 0, c);
  acb_mul_si(c, A22(Aij), -2, prec);
  acb_poly_set_coeff_acb(P1, 1, c);  

  acb_poly_set_coeff_acb(P2, 0, A11(Aij));
  acb_neg(c, A22(Aij));
  acb_poly_set_coeff_acb(P2, 2, c);

  acb_poly_set_coeff_acb(P3, 0, A11(Aij));
  acb_mul_si(c, A12(Aij), 2, prec);
  acb_poly_set_coeff_acb(P3, 1, c);
  acb_poly_set_coeff_acb(P3, 2, A22(Aij));

  /* Set crv */
  acb_poly_zero(crv);
  
  acb_poly_pow_ui(term, P1, 3, prec);
  acb_mul(c, A33(Aij), a111(aijk), prec);
  acb_neg(c, c);
  acb_poly_scalar_mul(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);

  acb_poly_pow_ui(term, P1, 2, prec);
  acb_poly_mul(term, term, P2, prec);
  acb_mul(c, A33(Aij), a112(aijk), prec);
  acb_mul_si(c, c, -3, prec);
  acb_poly_scalar_mul(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);

  acb_poly_pow_ui(term, P2, 2, prec);
  acb_poly_mul(term, term, P1, prec);
  acb_mul(c, A33(Aij), a122(aijk), prec);
  acb_mul_si(c, c, -3, prec);
  acb_poly_scalar_mul(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);
  
  acb_poly_pow_ui(term, P3, 2, prec);
  acb_poly_mul(term, term, P1, prec);
  acb_mul(c, A22(Aij), a133(aijk), prec);
  acb_mul_si(c, c, 3, prec);
  acb_poly_scalar_mul(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);
  
  acb_poly_pow_ui(term, P2, 3, prec);
  acb_mul(c, A33(Aij), a222(aijk), prec);
  acb_neg(c, c);
  acb_poly_scalar_mul(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);
  
  acb_poly_pow_ui(term, P3, 2, prec);
  acb_poly_mul(term, term, P2, prec);
  acb_mul(c, A22(Aij), a233(aijk), prec);
  acb_mul_si(c, c, 3, prec);
  acb_poly_scalar_mul(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);
    
  _acb_vec_clear(ABCD, 4);
  _acb_vec_clear(Aij, 4);
  _acb_vec_clear(aijk, 6);
  acb_poly_clear(P1);
  acb_poly_clear(P2);
  acb_poly_clear(P3);
  acb_clear(c);
  acb_poly_clear(term);
}

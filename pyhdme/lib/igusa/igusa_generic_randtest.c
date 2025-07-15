
#include "igusa.h"

void igusa_generic_randtest(acb_poly_t crv, acb_ptr I, flint_rand_t state, slong prec)
{
  acb_ptr roots;
  int stop = 0;

  roots = _acb_vec_init(6);
  while (!stop)
    {
      acb_randtest_precise(&roots[0], state, prec, 1);
      acb_randtest_precise(&roots[1], state, prec, 1);
      acb_randtest_precise(&roots[2], state, prec, 1);
      acb_randtest_precise(&roots[3], state, prec, 1);
      acb_randtest_precise(&roots[4], state, prec, 1);
      acb_randtest_precise(&roots[5], state, prec, 1);
      
      acb_poly_product_roots(crv, roots, 6, prec);
      igusa_from_curve(I, crv, prec);
      igusa_IC(I, I, prec);
      if (igusa_has_generic_automorphisms(I, prec)) stop = 1;
    }

  _acb_vec_clear(roots, 6);
}

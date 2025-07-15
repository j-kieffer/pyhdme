
#include "igusa.h"

void thomae_rosenhain(acb_ptr ros, acb_srcptr roots, slong prec)
{
  acb_ptr proj_roots;
  acb_t den;
  slong k;
  
  proj_roots = _acb_vec_init(5);
  acb_init(den);
  
  /* Bring root 6 to infinity */
  for (k = 0; k < 5; k++)
    {
      acb_sub(&proj_roots[k], &roots[k], &roots[5], prec);
      acb_inv(&proj_roots[k], &proj_roots[k], prec);
    }
  /* Bring roots 1 and 2 to 0, 1 */
  acb_sub(den, &proj_roots[1], &proj_roots[0], prec);
  acb_inv(den, den, prec);
  for (k = 2; k < 5; k++)
    {
      acb_sub(&proj_roots[k], &proj_roots[k], &proj_roots[0], prec);
      acb_mul(&proj_roots[k], &proj_roots[k], den, prec);
    }
  /* Set Rosenhain invariants: l, m, n */
  acb_set(&ros[0], &proj_roots[2]);
  acb_set(&ros[1], &proj_roots[3]);
  acb_set(&ros[2], &proj_roots[4]); 

  _acb_vec_clear(proj_roots, 5);
  acb_clear(den);
}


#include "siegel.h"

static int
not_minkowski_reduced_g2(const arb_mat_t r, slong prec)
{
  arb_t lhs, rhs;
  int res = 0;
  
  arb_init(lhs);
  arb_init(rhs);
  
  /* NOT 0 \leq 2y_3 */
  arb_mul_si(rhs, arb_mat_entry(r, 0, 1), 2, prec);
  res = arb_is_negative(rhs);
  
  /* NOT 2y_3 \leq y_1 */
  arb_mul_si(lhs, arb_mat_entry(r, 0, 1), 2, prec);
  arb_set(rhs, arb_mat_entry(r, 0, 0));
  res = res || arb_gt(lhs, rhs);

  /* NOT y_1 \leq y_2 */
  arb_set(lhs, arb_mat_entry(r, 0, 0));
  arb_set(rhs, arb_mat_entry(r, 1, 1));
  res = res || arb_gt(lhs, rhs);

  /* NOT y1 > 0 */
  res = res || arb_is_nonpositive(arb_mat_entry(r, 0, 0));
  /* NOT y2 > 0 */
  res = res || arb_is_nonpositive(arb_mat_entry(r, 1, 1));
  /* NOT det > 0 */
  arb_mat_det(lhs, r, prec);
  res = res || arb_is_nonpositive(lhs);

  arb_clear(lhs);
  arb_clear(rhs);
  return res;
}

int
arb_mat_not_minkowski_reduced(const arb_mat_t r, slong prec)
{
  slong g = arb_mat_nrows(r);
  int res;
  switch (g)
    {
    case 1:
      res = arb_is_nonpositive(arb_mat_entry(r, 0, 0));
      break;
    case 2:
      res = not_minkowski_reduced_g2(r, prec);
      break;
    default:
      flint_fprintf(stderr, "Error: Minkowski reduction test not implemented for g=%wd\n", g);
      flint_abort();
    }
  return res;
}

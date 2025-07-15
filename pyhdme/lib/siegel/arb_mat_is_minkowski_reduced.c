
#include "siegel.h"

static int
is_minkowski_reduced_g2(const arb_mat_t r, const arb_t tol, slong prec)
{
  arb_t lhs, rhs;
  int res = 1;
  
  arb_init(lhs);
  arb_init(rhs);
  
  /* -ey_1 \leq 2y_3 */
  arb_mul(lhs, tol, arb_mat_entry(r, 0, 0), prec);
  arb_neg(lhs, lhs);
  arb_mul_si(rhs, arb_mat_entry(r, 0, 1), 2, prec);
  res = arb_lt(lhs, rhs);
  
  /* 2y_3 \leq (1+e)y_1 */
  arb_mul_si(lhs, arb_mat_entry(r, 0, 1), 2, prec);
  arb_add_si(rhs, tol, 1, prec);
  /*arb_printd(rhs, 300), flint_printf("\n");*/
  arb_mul(rhs, rhs, arb_mat_entry(r, 0, 0), prec);
  res = res && arb_lt(lhs, rhs);

  arb_sub(lhs, rhs, lhs, prec);
  /*arb_printd(lhs, 30);flint_printf("\n");
    arb_printd(tol, 30);flint_printf("\n");
    flint_printf("res=%d\n", res);*/

  /* y_1 \leq (1+e)y_3 */
  arb_set(lhs, arb_mat_entry(r, 0, 0));
  arb_add_si(rhs, tol, 1, prec);
  arb_mul(rhs, rhs, arb_mat_entry(r, 1, 1), prec);
  res = res && arb_lt(lhs, rhs);

  arb_clear(lhs);
  arb_clear(rhs);
  return res;
}

int
arb_mat_is_minkowski_reduced(const arb_mat_t r, const arb_t tol, slong prec)
{
  slong g = arb_mat_nrows(r);
  int res;
  switch (g)
    {
    case 1:
      res = arb_is_positive(arb_mat_entry(r, 0, 0));
      break;
    case 2:
      res = is_minkowski_reduced_g2(r, tol, prec);
      break;
    default:
      flint_fprintf(stderr, "Error: Minkowski reduction test not implemented for g=%wd\n", g);
      flint_abort();
    }
  return res;
}

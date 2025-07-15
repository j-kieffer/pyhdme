
#include "siegel.h"

static void
arb_mat_lambda_g2(arb_t lambda, const arb_mat_t m, slong prec)
{
  arb_t tr;
  arb_t det;
  
  arb_init(tr);
  arb_init(det);

  arb_mat_trace(tr, m, prec);
  arb_mat_det(det, m, prec);

  arb_sqr(lambda, tr, prec);
  arb_submul_si(lambda, det, 4, prec);
  arb_nonnegative_part(lambda, lambda);
  
  arb_sqrt(lambda, lambda, prec);
  arb_sub(lambda, tr, lambda, prec);
  arb_mul_2exp_si(lambda, lambda, -1);

  arb_clear(tr);
  arb_clear(det);
}

static void
arb_mat_lambda_gg(arb_t lambda, const arb_mat_t m, slong prec)
{
  flint_fprintf(stderr, "\n(arb_mat_lambda_gg) Not implemented\n");
  flint_abort();
}

void
arb_mat_lambda(arb_t lambda, const arb_mat_t m, slong prec)
{
  slong g = arb_mat_nrows(m);
  switch(g)
    {
    case 1:
      arb_set(lambda, arb_mat_entry(m, 0, 0));
      break;
    case 2:
      arb_mat_lambda_g2(lambda, m, prec);
      break;
    default:
      arb_mat_lambda_gg(lambda, m, prec);
      break;
    }
}

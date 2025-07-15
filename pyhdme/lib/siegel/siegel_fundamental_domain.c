
#include "siegel.h"
#include "verbose.h"

int siegel_fundamental_domain(acb_mat_t w, fmpz_mat_t m,
    const acb_mat_t z, const arb_t tol, slong prec)
{
  int stop = 0;
  int res = 1;
  int v = get_siegel_verbose();
  int j;

  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t test_matrix;
  fmpz_mat_t step;
  fmpz_mat_t u;
  arb_mat_t im;
  acb_mat_t x;
  acb_mat_t star;
  acb_t det;
  arb_t absdet;

  fmpz_mat_init(test_matrix, 2*g, 2*g);
  fmpz_mat_init(step, 2*g, 2*g);
  fmpz_mat_init(u, g, g);
  arb_mat_init(im, g, g);
  acb_mat_init(x, g, g);
  acb_mat_init(star, g, g);
  acb_init(det);
  arb_init(absdet);

  fmpz_mat_one(m);

  while (!stop) {
    /* Update x */
    res = res && siegel_transform(x, m, z, prec);
    if (!res) {
      if (v) {
        flint_printf("\n(siegel_fundamental_domain) Precision lost in update:\n");
        flint_printf("m = "); fmpz_mat_print(m); flint_printf("\n\n");
        flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n\n");
        flint_printf("tol = "); arb_printd(tol, 5); flint_printf("\n");
      }
      break;
    }

    /* Reduce imaginary part */
    acb_mat_get_imag(im, x);
    res = arb_mat_minkowski_reduce(im, u, im, tol, prec);
    fmpz_mat_diagonal_symplectic(step, u);

    fmpz_mat_mul(m, step, m);
    res = res && siegel_transform(x, m, z, prec);
    if (!res) {
      if (v) {
      flint_printf("\n(siegel_fundamental_domain) Precision lost in Minkowski reduction:\n");
      flint_printf("im = "); arb_mat_printd(im, 30); flint_printf("\n\n");
      flint_printf("u = "); fmpz_mat_print_pretty(u); flint_printf("\n\n");
      flint_printf("tol = "); arb_printd(tol, 5); flint_printf("\n");
      }
      break;
    }

    /* Reduce real part */
    res = siegel_reduce_real(x, step, x, tol, prec);
    fmpz_mat_mul(m, step, m);
    if (!res) {
      if (v) {
        flint_printf("\n(siegel_fundamental_domain) Precision lost in real reduction:\n");
        flint_printf("x = "); acb_mat_printd(x, 30); flint_printf("\n\n");
        flint_printf("tol = "); arb_printd(tol, 5); flint_printf("\n");
      }
      break;
    }

    /* Apply test matrix */
    /* To be improved: compute determinants at smaller precision,
       select the smallest one */
    stop = 1;
    for (j = 0; j < siegel_nb_test_matrices(g); j++)
    {
      siegel_test_matrix(step, j);
      siegel_star(star, step, x, prec);
      acb_mat_det(det, star, prec);
      acb_abs(absdet, det, prec);
      arb_sub_si(absdet, absdet, 1, prec);
      if (arb_is_negative(absdet))
      {
        stop = 0;
        fmpz_mat_mul(m, step, m);
        break; /* for loop */
      }
    }
  }

  res = res && siegel_transform(w, m, z, prec);

  fmpz_mat_clear(test_matrix);
  fmpz_mat_clear(step);
  fmpz_mat_clear(u);
  arb_mat_clear(im);
  acb_mat_clear(x);
  acb_mat_clear(star);
  acb_clear(det);
  arb_clear(absdet);
  return res;
}


#include "hecke.h"

int hecke_charpoly(fmpz_poly_t pol, slong ell, slong wt)
{
  slong weights[4] = IGUSA_WEIGHTS;
  slong nb = cov_nb_monomials(wt, 4, weights);
  acb_mat_struct* pts;
  acb_mat_t basis_inv;
  acb_mat_t image;
  acb_poly_t charpoly;
  arf_t rad;
  slong gap;
  slong prec = hecke_charpoly_startprec(wt);
  int v = get_hecke_verbose();
  slong k;
  int res = 0;

  pts = flint_malloc(nb * sizeof(acb_mat_struct));
  for (k = 0; k < nb; k++) acb_mat_init(&pts[k], 2, 2);
  acb_mat_init(basis_inv, nb, nb);
  acb_mat_init(image, nb, nb);
  acb_poly_init(charpoly);
  arf_init(rad);

  while (!res) {
    if (v) flint_printf("(hecke_charpoly) Start new run at precision %wd\n", prec);

    res = hecke_select_basepts(pts, basis_inv, wt, prec);
    if (res) res = hecke_image_matrix(image, nb, pts, wt, ell, prec);
    if (res) {
      acb_mat_mul(image, image, basis_inv, prec);
      acb_mat_charpoly(charpoly, image, prec);
      res = acb_poly_round(pol, rad, charpoly, nb);
      gap = arf_abs_bound_lt_2exp_si(rad);
      prec += gap/2 + 50;
    }

    prec = hecke_charpoly_nextprec(prec);

    if (prec > n_pow(10, 6))
    {
      flint_printf("(hecke_charpoly) Precision too high, abort.\n");
      fflush(stdout);
      flint_abort();
    }
  }

  for (k = 0; k < nb; k++) acb_mat_clear(&pts[k]);
  flint_free(pts);
  acb_mat_clear(basis_inv);
  acb_mat_clear(image);
  acb_poly_clear(charpoly);
  arf_clear(rad);
  return res;
}

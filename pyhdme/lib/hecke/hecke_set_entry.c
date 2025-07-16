
#include <flint/acb_theta.h>
#include "hecke.h"

int hecke_set_entry(hecke_t H, slong k, const fmpz_mat_t gamma, slong prec)
{
  acb_mat_t im;
  fmpz_mat_t eta;
  int res;
  slong nb = hecke_nb(H);

  acb_mat_init(im, 2, 2);
  fmpz_mat_init(eta, 4, 4);

  /* Sanity check before using macros */
  if (k < 0 || k >= nb) {
    flint_printf("(hecke_set_entry) Error: illegal index %wd, expected between 0 and %wd\n",
        k, nb-1);
    fflush(stdout);
    flint_abort();
  }

  acb_siegel_transform(im, gamma, hecke_tau(H), prec);
  acb_siegel_reduce(eta, im, prec);
  acb_siegel_transform(hecke_isog(H, k), eta, im, prec);
  res = acb_siegel_is_reduced(hecke_isog(H, k), -HECKE_RED_TOL_BITS, prec);

  if (res) {
    /* Can fill hecke_coset, hecke_star */
    fmpz_mat_mul(hecke_coset(H, k), eta, gamma);
    acb_siegel_cocycle(hecke_star(H, k), hecke_coset(H, k), hecke_tau(H), prec);
    acb_mat_det(hecke_stardet(H, k), hecke_star(H, k), prec);

    /* Compute projective vector of theta constants */
    res = theta2_unif(hecke_theta2(H, k), hecke_isog(H, k), prec);
  }
  if (res) res = theta2_renormalize(hecke_theta2(H, k), hecke_theta2(H, k), prec);
  if (res) igusa_from_theta2(hecke_I(H, k), hecke_theta2(H, k), prec);

  acb_mat_clear(im);
  fmpz_mat_clear(eta);
  return res;
}

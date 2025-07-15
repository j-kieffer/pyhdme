
#include "hecke.h"

int hecke_image_matrix(acb_mat_t image, slong nb, const acb_mat_struct* pts,
    slong wt, slong ell, slong prec)
{
  hecke_t H;
  acb_ptr act;
  slong weights[4] = IGUSA_WEIGHTS;

  slong m = siegel_nb_cosets(ell);
  slong n = cov_nb_monomials(wt, 4, weights);
  int res = 1;

  hecke_init(H, m);
  act = _acb_vec_init(n);

  for (slong k = 0; k < nb; ++k) {
    res = hecke_set_tau(H, &pts[k], prec);
    if (!res) break;
    res = hecke_collect_siegel(H, ell, prec);
    if (!res) break;
    hecke_action_all_monomials(act, H, wt, ell, prec);
    for (slong j = 0; j < n; ++j) {
      acb_set(acb_mat_entry(image, j, k), &act[j]);
    }
  }

  hecke_clear(H);
  _acb_vec_clear(act, n);
  return res;
}


#include "hecke.h"

/* See also hecke_collect_T1 */

int hecke_collect_siegel(hecke_t H, slong ell, slong prec) {
  time_pair start; timestamp_mark(&start);
  slong nb = hecke_nb(H);
  int res = 1;
  int v = get_hecke_verbose();
  slong steps_completed = 0;

  hecke_ell(H) = ell;
  hecke_check_nb(H, siegel_nb_cosets(ell));

  if (v) hecke_collect_verbose_start(nb);

#pragma omp parallel
  {
    int master = 0;
#pragma omp single
    {
      master = 1;
    }
    fmpz_mat_t gamma;
    fmpz_mat_init(gamma, 4, 4);
#pragma omp for schedule(static)
    /* Loop over all cosets to compute desired data */
    for (slong k = 0; k < nb; ++k) {
      if (!res) continue; // OpenMP does not allow break, and thus we continue
      if (master && v) progress_bar(steps_completed, nb, "(hecke_collect)");

      siegel_coset(gamma, k, ell);
      int localres = hecke_set_entry(H, k, gamma, prec);
#pragma omp atomic
      res &= localres;
#pragma omp atomic
      ++steps_completed;
    }

    fmpz_mat_clear(gamma);
    if (master && v) {
      progress_bar(steps_completed, nb, "(hecke_collect)");
      if (!res)
        flint_printf("(hecke_collect) Warning: computation aborted due to low precision\n");
    }
  }

  /* Set normalization factor */
  hecke_norm_ind(H) = n_pow(ell, 2);
  fmpz_set_si(hecke_norm_all(H), ell);
  fmpz_pow_ui(hecke_norm_all(H), hecke_norm_all(H), 2*hecke_nb(H) - (ell*ell + ell + 2));
  hecke_prod_ec(H) = n_pow(ell+1, 2);

  if (v) report_end(start);
  return res;
}

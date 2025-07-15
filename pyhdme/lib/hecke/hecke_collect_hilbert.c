
#include "hecke.h"

int hecke_collect_hilbert(hecke_t H, const fmpz_poly_t beta,
			  slong ell, slong delta, slong prec)
{  
  slong k;
  fmpz_poly_mat_t m;
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  int res = 1;
  int v = get_hecke_verbose();

  fmpz_poly_mat_init(m, 2, 2);
  fmpz_mat_init(gamma, 4, 4);

  /* Set Hecke context */
  hecke_ell(H) = ell;
  fmpz_poly_set(hecke_beta(H), beta);
  hecke_check_nb(H, hilbert_nb_cosets(ell, delta));
  
  if (v) hecke_collect_verbose_start(nb);

  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      if (v) hecke_collect_print_status(res, k, nb);
      if (!res) break;
      
      hilbert_coset(m, k, beta, ell, delta);
      hilbert_mat_map(gamma, m, delta);
      fmpz_mat_mul(gamma, gamma, hecke_eta(H));

      res = hecke_set_entry(H, k, gamma, prec);  
    }
  
  if (v) flint_printf("\n");

  hecke_norm_ind(H) = ell;
  fmpz_set_si(hecke_norm_all(H), ell);
  fmpz_pow_ui(hecke_norm_all(H), hecke_norm_all(H), hecke_nb(H) - 1);
  hecke_prod_ec(H) = 0;

  fmpz_poly_mat_clear(m);
  fmpz_mat_clear(gamma);
  return res;  
}


#include "hecke.h"

int hecke_collect_hilbert_sym(hecke_t H,
			      slong ell, slong delta, slong prec)
{
  slong k;
  fmpz_poly_mat_t m;
  fmpz_mat_t gamma;
  fmpz_poly_t beta, betabar;
  slong nb = hecke_nb(H);
  int res = 1;
  int v = get_hecke_verbose();

  fmpz_poly_mat_init(m, 2, 2);
  fmpz_mat_init(gamma, 4, 4);
  fmpz_poly_init(beta);
  fmpz_poly_init(betabar);

  /* Set Hecke context */
  hecke_ell(H) = ell;
  res = hilbert_splits(beta, ell, delta);
  if (!res)
    {
      flint_printf("(hecke_collect_hilbert_sym) %wd does not split", ell);
      fflush(stdout);
      flint_abort();
    }
  
  fmpz_poly_set(hecke_beta(H), beta);
  hilbert_conjugate(betabar, beta, delta);
  hecke_check_nb(H, 2*hilbert_nb_cosets(ell, delta));
  
  if (v) hecke_collect_verbose_start(nb);
  
  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      if (v) hecke_collect_print_status(res, k, nb);
      if (!res) break;

      if (k < hilbert_nb_cosets(ell, delta))
	{
	  hilbert_coset(m, k, beta, ell, delta);
	}
      else
	{
	  hilbert_coset(m, k - hilbert_nb_cosets(ell, delta),
			betabar, ell, delta);
	}
      
      hilbert_mat_map(gamma, m, delta);
      fmpz_mat_mul(gamma, gamma, hecke_eta(H));
      
      res = hecke_set_entry(H, k, gamma, prec);    
    }
  
  if (v) flint_printf("\n");

  hecke_norm_ind(H) = ell;
  fmpz_set_si(hecke_norm_all(H), ell);
  fmpz_pow_ui(hecke_norm_all(H), hecke_norm_all(H), hecke_nb(H) - 2);
  hecke_prod_ec(H) = 0;

  fmpz_poly_mat_clear(m);
  fmpz_mat_clear(gamma);
  fmpz_poly_clear(beta);
  fmpz_poly_clear(betabar);
  return res;  
}


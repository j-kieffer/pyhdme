
#include "hecke.h"

int hecke_collect_T1_with_line(hecke_t H, const fmpz_mat_t L,
			       slong ell, slong prec)
{
  slong k, j;
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  slong nb_all = siegel_nb_T1_cosets(ell);
  
  int res = 1;
  int v = get_hecke_verbose();
  
  fmpz_mat_init(gamma, 4, 4);  
  hecke_ell(H) = ell;
  hecke_check_nb(H, siegel_nb_T1_cosets_with_line(ell));
    
  if (v) hecke_collect_verbose_start(nb);

  /* Loop over all cosets to compute desired data */
  k = 0;
  for (j = 0; j < nb_all; j++)
    {      
      siegel_T1_coset(gamma, j, ell);
      if (siegel_T1_coset_contains_line_dual(gamma, L, ell))
	{	  
	  if (v) hecke_collect_print_status(res, k, nb);
	  if (!res) break;
	  
	  res = hecke_set_entry(H, k, gamma, prec);
	  k++;
	}
    }
  if (v) flint_printf("\n");

  hecke_norm_ind(H) = n_pow(ell, 3);
  fmpz_set_si(hecke_norm_all(H), ell);
  fmpz_pow_ui(hecke_norm_all(H), hecke_norm_all(H), 3*hecke_nb(H));
  hecke_prod_ec(H) = 0;
  
  fmpz_mat_clear(gamma);
  return res;  
}

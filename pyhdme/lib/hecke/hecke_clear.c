
#include "hecke.h"

void hecke_clear(hecke_t H)
{
  slong k;
  slong nb = hecke_nb(H);
  
  acb_mat_clear(hecke_tau(H));
  _acb_vec_clear(hecke_theta2_tau(H), 16);
  _acb_vec_clear(hecke_I_tau(H), 4);
  _acb_vec_clear(hecke_t1t2(H), 2);
  fmpz_mat_clear(hecke_eta(H));
  fmpz_poly_clear(hecke_beta(H));
  fmpz_clear(hecke_norm_all(H));
  
  for (k = 0; k < nb; k++) fmpz_mat_clear(hecke_coset(H, k));
  for (k = 0; k < nb; k++) acb_mat_clear(hecke_star(H, k));
  for (k = 0; k < nb; k++) acb_mat_clear(hecke_isog(H, k));  

  flint_free(H->cosets);
  flint_free(H->isog);
  flint_free(H->stars);
  _acb_vec_clear(H->stardets, nb);
  _acb_vec_clear(H->theta2, 16*nb);
  _acb_vec_clear(H->I, 4*nb);
  _acb_vec_clear(H->I_norm, 4*nb);
}

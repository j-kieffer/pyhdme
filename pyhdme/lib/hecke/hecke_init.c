
#include "hecke.h"

void hecke_init(hecke_t H, slong nb)
{
  slong k;
  
  acb_mat_init(hecke_tau(H), 2, 2);
  hecke_theta2_tau(H) = _acb_vec_init(16);
  hecke_I_tau(H) = _acb_vec_init(4);
  hecke_t1t2(H) = _acb_vec_init(2);
  fmpz_mat_init(hecke_eta(H), 4, 4);
  fmpz_poly_init(hecke_beta(H));

  hecke_nb(H) = nb;
  fmpz_init(hecke_norm_all(H));
  H->cosets = flint_malloc(nb * sizeof(fmpz_mat_struct));  
  H->isog = flint_malloc(nb * sizeof(acb_mat_struct));
  H->stars = flint_malloc(nb * sizeof(acb_mat_struct));
  H->stardets = _acb_vec_init(nb);
  H->theta2 = _acb_vec_init(16*nb);
  H->I = _acb_vec_init(4*nb);
  H->I_norm = _acb_vec_init(4*nb);
  
  for (k = 0; k < nb; k++) fmpz_mat_init(hecke_coset(H, k), 4, 4);
  for (k = 0; k < nb; k++) acb_mat_init(hecke_star(H, k), 2, 2);
  for (k = 0; k < nb; k++) acb_mat_init(hecke_isog(H, k), 2, 2);  
}

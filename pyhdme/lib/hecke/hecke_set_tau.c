
#include "hecke.h"

int hecke_set_tau(hecke_t H, const acb_mat_t tau, slong prec)
{
  int res;
  int v = get_hecke_verbose();

  hecke_prec(H) = prec;
  acb_mat_set(hecke_tau(H), tau);
  
  res = theta2_unif(hecke_theta2_tau(H), tau, prec);
  if (res) theta2_renormalize(hecke_theta2_tau(H), hecke_theta2_tau(H), prec);
  if (res) igusa_from_theta2(hecke_I_tau(H), hecke_theta2_tau(H), prec);
  if (v && !res)
    {
      flint_printf("(hecke_set_tau) Warning: computation aborted due to low precision\n");
    }
  
  return res;
}

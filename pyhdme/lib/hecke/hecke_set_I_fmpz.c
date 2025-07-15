
#include "hecke.h"

int hecke_set_I_fmpz(hecke_t H, fmpz* I, slong prec)
{
  int res;
  int v = get_hecke_verbose();

  hecke_prec(H) = prec;

  if (v) flint_printf("(hecke_set_I_fmpz) Computing period matrix...\n");
  res = tau_theta2_from_igusa_fmpz(hecke_tau(H), hecke_theta2_tau(H),
				   I, prec);
  if (v && res) flint_printf("(hecke_set_I_fmpz) Period matrix found.\n");
  
  if (res) igusa_from_theta2(hecke_I_tau(H), hecke_theta2_tau(H), prec);
  if (v && !res)
    {
      flint_printf("(hecke_set_I_fmpz) Warning: computation aborted due to low precision\n");
    }  
  
  return res;
}

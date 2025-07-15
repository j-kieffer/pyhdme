
#include "igusa.h"

int tau_theta2_from_igusa_ec(acb_mat_t tau, acb_ptr th2, fmpz* I, slong prec)
{
  acb_ptr j;
  slong k;
  int res = 1;

  j = _acb_vec_init(2);
  
  igusa_ec_j1j2(j, I, prec);  
  acb_mat_zero(tau);
  
  for (k = 0; k < 2; k++) res = res && igusa_ec_period(&j[k], &j[k], prec);
  if (arb_lt(acb_imagref(&j[1]), acb_imagref(&j[0]))) acb_swap(&j[0], &j[1]);
  for (k = 0; k < 2; k++) acb_set(acb_mat_entry(tau, k, k), &j[k]);
  if (res) res = theta2_unif(th2, tau, prec);
  if (res) res = theta2_renormalize(th2, th2, prec);

  _acb_vec_clear(j, 2);
  return res;
}

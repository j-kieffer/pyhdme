
#include "hecke.h"

/* See also hecke_set_I_fmpz */

int hecke_set_I_fmpz_with_lowprec(hecke_t H, fmpz* I, acb_srcptr th2_lp, slong prec) {
  int res;
  int v = get_hecke_verbose();
  time_pair start; timestamp_mark(&start);

  hecke_prec(H) = prec;

  if (v) flint_printf("(hecke_set_I_fmpz_with_lowprec) Refining period matrix...\n");
  res = tau_theta2_from_igusa_fmpz_with_lowprec(
      hecke_tau(H),
      hecke_theta2_tau(H),
      I,
      th2_lp,
      prec);
  if (v && res) flint_printf("(hecke_set_I_fmpz_with_lowprec) Period matrix found.\n");

  if (res) igusa_from_theta2(hecke_I_tau(H), hecke_theta2_tau(H), prec);
  if (v && !res) {
    flint_printf("(hecke_set_I_fmpz_with_lowprec) Warning: computation aborted due to low precision\n");
  }

  if (v) report_end(start);
  return res;
}

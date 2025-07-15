
#include "modular.h"

int siegel_direct_isog_Q(slong* nb, fmpz* all_I, fmpz* I, slong ell) {
  time_pair start; timestamp_mark(&start);

  slong prec = siegel_direct_isog_startprec(I, ell);
  hecke_t H;
  int stop = 0;
  int res = 1;
  int v = get_modeq_verbose();
  *nb = 0; /* we might never call hecke_all_isog_Q */

  hecke_init(H, siegel_nb_cosets(ell));

  while (!stop) {
    if (v) flint_printf("(siegel_direct_isog_Q) trying with prec = %wd\n", prec);
    res = hecke_set_I_fmpz(H, I, prec);

    if (res) res = hecke_collect_siegel(H, ell, prec);
    if (res) hecke_make_integral(H, I, prec);
    if (res) res = hecke_has_integral_precision(H, prec);
    if (res) res = hecke_all_isog_Q(nb, all_I, H, I, prec);

    prec = modeq_nextprec_generic(prec);
    stop = modeq_stop(res, prec);
  }
  if (v) flint_printf("(siegel_direct_isog_Q) succeeded with prec = %wd\n", prec);
  hecke_clear(H);
  if (get_hecke_verbose()) report_end(start);
  return res;
}

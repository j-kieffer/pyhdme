
#include "modular.h"

int siegel_2step_direct_isog_Q(slong* nb, fmpz* all_I, fmpz* I, slong ell) {
  time_pair start; timestamp_mark(&start);

  slong prec = siegel_direct_isog_startprec(I, ell);
  hecke_t H;
  int stop = 0;
  int res = 1;
  int v = get_modeq_verbose();
  *nb = 0; /* we might never call hecke_all_isog_Q */

  if (v) {
    flint_printf("(siegel_2step_direct_isog_Q) ell = %wd\n", ell);
    flint_printf("(siegel_2step_direct_isog_Q) I = ");
    _fmpz_vec_print(I, 4);
    flint_printf("\n");
  }
  hecke_init(H, siegel_nb_T1_cosets(ell));

  while (!stop) {
    if (v) flint_printf("(siegel_2step_direct_isog_Q) trying with prec = %wd\n", prec);
    res = hecke_set_I_fmpz(H, I, prec);

    if (res) res = hecke_collect_T1(H, ell, prec);
    if (res) hecke_make_integral(H, I, prec);
    if (res) res = hecke_has_integral_precision(H, prec);
    if (res) res = hecke_all_isog_Q(nb, all_I, H, I, prec);

    prec = modeq_nextprec_generic(prec);
    stop = modeq_stop(res, prec);
  }
  if (v) {
    flint_printf("(siegel_2step_direct_isog_Q) succeeded with prec = %wd\n", prec);
    flint_printf("(siegel_2step_direct_isog_Q) nb = %wd\n", *nb);
    if (*nb > 0) {
      flint_printf("(siegel_2step_direct_isog_Q) roots = [\n");
      for(slong i=0; i < *nb; ++i) {
        _fmpz_vec_print(&all_I[i*4], 4);
        flint_printf("\n");
      }
      flint_printf("]\n");
    } else {
      flint_printf("(siegel_2step_direct_isog_Q) roots = []\n");
    }
  }



  hecke_clear(H);
  if (v) report_end(start);
  return res;
}

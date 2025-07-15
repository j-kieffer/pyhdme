
#include "theta.h"

slong theta_newton_k2(acb_mat_t w, const acb_mat_t z, slong prec)
{
  arb_t y1;
  arb_t y2;
  arb_t t;
  arf_t u;
  fmpz_t v;
  slong k2;

  arb_init(y1);
  arb_init(y2);
  arb_init(t);
  arf_init(u);
  fmpz_init(v);
  
  arb_set(y1, acb_imagref(acb_mat_entry(z, 0, 0)));
  arb_set(y2, acb_imagref(acb_mat_entry(z, 1, 1)));

  arb_set_si(t, prec);
  arb_div_si(t, t, THETA_NEWTON_Y1, prec);
  arb_div(t, t, y1, prec);
  arb_log_base_ui(t, t, 2, prec);
  arb_get_ubound_arf(u, t, prec);
  arf_get_fmpz(v, u, ARF_RND_CEIL);
  k2 = fmpz_get_si(v);

  arb_div(t, y2, y1, prec);
  arb_log_base_ui(t, t, 4, prec);
  arb_get_ubound_arf(u, t, prec);
  arf_get_fmpz(v, u, ARF_RND_CEIL);
  k2 = FLINT_MIN(k2, fmpz_get_si(v) - 1);
  k2 = FLINT_MAX(k2, 0);

  acb_mat_set(w, z);
  acb_mul_2exp_si(acb_mat_entry(w, 0, 0), acb_mat_entry(w, 0, 0), k2);
  acb_mul_2exp_si(acb_mat_entry(w, 1, 1), acb_mat_entry(w, 1, 1), -k2);

  arb_clear(y1);
  arb_clear(y2);
  arb_clear(t);
  arf_clear(u);
  fmpz_clear(v);
  return k2;
}

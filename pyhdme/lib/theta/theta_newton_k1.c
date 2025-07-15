
#include "theta.h"

slong theta_newton_k1(acb_mat_t w, const acb_mat_t z, slong prec)
{
  arb_t abs_z1;
  arb_t abs_z2;
  arb_t t;
  arf_t u;
  fmpz_t v;
  slong k1;

  arb_init(abs_z1);
  arb_init(abs_z2);
  arb_init(t);
  arf_init(u);
  fmpz_init(v);

  acb_abs(abs_z1, acb_mat_entry(z, 0, 0), prec);
  acb_abs(abs_z2, acb_mat_entry(z, 1, 1), prec);
  arb_min(t, abs_z1, abs_z2, prec);
  arb_log_base_ui(t, t, 2, prec);
  arb_get_lbound_arf(u, t, prec);
  arf_get_fmpz(v, u, ARF_RND_FLOOR);
  k1 = fmpz_get_si(v);
  k1 = FLINT_MAX(k1, 0);

  /* If y1 is large enough, set k1=0 */
  arb_set(t, acb_imagref(acb_mat_entry(z, 0, 0)));
  arb_mul_si(t, t, THETA_NEWTON_Y1, prec);
  arb_sub_si(t, t, prec, prec);
  if (arb_is_positive(t)) k1 = 0;

  acb_mat_scalar_mul_2exp_si(w, z, -k1);

  arb_clear(abs_z1);
  arb_clear(abs_z2);
  arb_clear(t);
  arf_clear(u);
  fmpz_clear(v);

  return k1;
}

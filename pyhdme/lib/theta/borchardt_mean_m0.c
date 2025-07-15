
#include "theta.h"

/* Since 1 is one of the initial values, we can choose the complex
   argument with a discontinuity along ]-infty,0] */
/* This computation may not be optimal. */

void
borchardt_mean_m0(arb_t m0, acb_srcptr a, slong prec)
{
  arb_ptr args;
  acb_ptr a_scaled;
  arf_t tmax, tmin, t;
  acb_t z;
  int i;

  args = _arb_vec_init(4);
  a_scaled = _acb_vec_init(4);
  arf_init(tmax);
  arf_init(tmin);
  arf_init(t);
  acb_init(z);

  /* Compute arguments */
  for (i = 0; i < 4; i++) acb_arg(&args[i], &a[i], prec);

  /* Compute bounds tmin, tmax */
  arb_get_lbound_arf(tmin, &args[0], prec);
  arb_get_ubound_arf(tmax, &args[0], prec);
  for (i = 1; i < 4; i++)
    {
      arb_get_lbound_arf(t, &args[i], prec);
      arf_min(tmin, tmin, t);
      arb_get_ubound_arf(t, &args[i], prec);
      arf_max(tmax, tmax, t);
    }

  /* Scale by z = exp(-it) where t is the midpoint */
  arf_add(t, tmin, tmax, prec, ARF_RND_NEAR); /* Rounding is unimportant here */
  arf_mul_2exp_si(t, t, -1);
  arb_set_arf(m0, t);
  acb_set_arb(z, m0);
  acb_neg(z, z);
  acb_mul_onei(z, z);
  acb_exp(z, z, prec);
  _acb_vec_scalar_mul(a_scaled, a, 4, z, prec);

  /* flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
     for (i = 0; i < 4; i++)
     {
     flint_printf("a[%wd] = ", i); acb_printd(&a[i], 30); flint_printf("\n");
     flint_printf("a_scaled[%wd] = ", i); acb_printd(&a_scaled[i], 30); flint_printf("\n");
     } */

  /* Set tmin to the minimum of real values */
  arb_get_lbound_arf(tmin, acb_realref(&a_scaled[0]), prec);
  for (i = 1; i < 4; i++)
    {
      arb_get_lbound_arf(t, acb_realref(&a_scaled[i]), prec);
      arf_min(tmin, tmin, t);
    }
  arb_set_arf(m0, tmin);

  _arb_vec_clear(args, 4);
  _acb_vec_clear(a_scaled, 4);
  arf_clear(tmax);
  arf_clear(tmin);
  arf_clear(t);
  acb_clear(z);
}

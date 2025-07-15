
#include "theta.h"

int
borchardt_mean_nb_steps_before_quad_conv(fmpz_t nb, acb_srcptr a, slong prec)
{
  arb_t m0, M0, Delta0, num, den;
  arf_t sup;
  slong g = 2;
  slong k;
  int res = 1;

  arb_init(m0);
  arb_init(M0);
  arb_init(Delta0);
  arb_init(num);
  arb_init(den);
  arf_init(sup);

  borchardt_mean_m0(m0, a, prec);
  borchardt_mean_M0(M0, a, prec);
  borchardt_mean_Delta0(Delta0, a, prec);

  if (!arb_is_positive(m0)
      || !arb_is_finite(m0)
      || !arb_is_finite(M0)
      || !arb_is_finite(Delta0)
      || arb_contains_negative(Delta0)) res = 0;

  if (!res && get_borchardt_verbose())
    {
      flint_printf("(borchardt_mean_nb_steps_before_quad_conv) Failure values:\n");
      flint_printf("m0 = "); arb_printd(m0, 30); flint_printf("\n");
      flint_printf("M0 = "); arb_printd(M0, 30); flint_printf("\n");
      flint_printf("Delta0 = "); arb_printd(Delta0, 30); flint_printf("\n");
      flint_printf("prec = %wd\n", prec);
      flint_printf("a: ");
      for (k = 0; k < 4; k++)
	{
	  acb_printd(&a[k], 30); flint_printf("\n");
	}
    }
  
  if (res)
    {      
      if (arb_is_zero(Delta0))
	{
	  fmpz_zero(nb);
	}
      else
	{
	  arb_div(num, m0, Delta0, prec);
	  arb_div_si(num, num, 7, prec);
	  arb_log(num, num, prec);
	  
	  arb_one(den);
	  arb_mul_2exp_si(den, den, -g);
	  arb_neg(den, den);
	  arb_add_si(den, den, 1, prec);
	  arb_log(den, den, prec);

	  arb_div(num, num, den, prec);
	  if (arb_is_negative(num)) arb_zero(num);
      
	  arb_get_ubound_arf(sup, num, prec);
	  arf_ceil(sup, sup);
	  arf_get_fmpz(nb, sup, ARF_RND_NEAR);
	}
    }
  
  arb_clear(m0);
  arb_clear(M0);
  arb_clear(Delta0);
  arb_clear(num);
  arb_clear(den);
  arf_clear(sup);
  return res;
}

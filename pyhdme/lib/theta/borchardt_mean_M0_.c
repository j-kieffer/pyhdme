
#include "theta.h"

void
borchardt_mean_M0(arb_t M0, acb_srcptr a, slong prec)
{
  arf_t abs;
  arf_t M0_arf;
  int i;

  arf_init(abs);
  arf_init(M0_arf);
  
  arf_zero(M0_arf);
  for (i = 0; i < 4; i++)
    {
      acb_get_abs_ubound_arf(abs, &a[i], prec);
      arf_max(M0_arf, M0_arf, abs);
    }
  arb_set_arf(M0, M0_arf);

  arf_clear(abs);
  arf_clear(M0_arf);
}

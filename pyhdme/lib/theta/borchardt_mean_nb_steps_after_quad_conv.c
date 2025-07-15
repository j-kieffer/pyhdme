
#include "theta.h"

void
borchardt_mean_nb_steps_after_quad_conv(fmpz_t nb, acb_srcptr a, slong prec)
{
  arb_t M0, num;
  arf_t sup;

  arb_init(M0);
  arb_init(num);
  arf_init(sup);

  borchardt_mean_M0(M0, a, prec);
  
  arb_div_si(num, M0, 7, prec);
  arb_log_base_ui(num, num, 2, prec);
  arb_add_si(num, num, prec, prec);
  arb_add_si(num, num, 1, prec);
  arb_log_base_ui(num, num, 2, prec);
  
  arb_get_ubound_arf(sup, num, prec);
  arf_ceil(sup, sup);
  arf_get_fmpz(nb, sup, ARF_RND_NEAR);
        
  arb_clear(M0);
  arb_clear(num);
  arf_clear(sup);
}

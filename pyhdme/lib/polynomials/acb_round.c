
#include "polynomials.h"

int acb_round(fmpz_t c, arf_t radius, const acb_t x)
{
  fmpz_t re, im;
  int res;
  slong radius_prec = HDME_RD_RADIUS_PREC;
  
  fmpz_init(re);
  fmpz_init(im);

  if (!acb_contains_int(x))
    {
      flint_printf("(acb_round) Error: contains no integer\n");
      acb_printd(x, 10); flint_printf("\n");
      fflush(stdout);
      flint_abort();
    }

  acb_get_rad_ubound_arf(radius, x, radius_prec);

  res = arb_get_unique_fmpz(re, acb_realref(x))
    && arb_get_unique_fmpz(im, acb_imagref(x));
  if (res) res = fmpz_is_zero(im);
  if (res) fmpz_set(c, re);

  fmpz_clear(re);
  fmpz_clear(im);
  return res;
}


#include "igusa.h"

int igusa_is_g2_curve_fmpz(fmpz* I)
{
  return !fmpz_is_zero(igusa_chi10(I));
}

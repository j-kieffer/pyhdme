
#include "igusa.h"

void igusa_streng(acb_ptr S, acb_srcptr I, slong prec)
{
  acb_mul_si(&S[0], igusa_psi4(I), 4, prec);
  acb_mul_si(&S[1], igusa_psi6(I), 4, prec);
  acb_mul_si(&S[2], igusa_chi10(I), -n_pow(2,12), prec);
  acb_mul_si(&S[3], igusa_chi12(I), n_pow(2,15), prec);
}

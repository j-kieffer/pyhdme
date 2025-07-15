
#include "igusa.h"

void igusa_streng_fmpz(fmpz* S, fmpz* I)
{
  slong weights[4] = IGUSA_HALFWEIGHTS;
  
  fmpz_mul_si(&S[0], &I[0], 4);
  fmpz_mul_si(&S[1], &I[1], 4);
  fmpz_mul_si(&S[2], &I[2], -n_pow(2,12));
  fmpz_mul_si(&S[3], &I[3], n_pow(2,15));
  
  cov_normalize_fmpz(S, S, 4, weights);
}

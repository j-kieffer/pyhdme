
#include "igusa.h"

void igusa_from_streng_fmpz(fmpz* I, fmpz* S)
{
  slong weights[4] = IGUSA_HALFWEIGHTS;
  
  cov_rescale_fmpz_si(I, S, 8, 4, weights);  
  fmpz_divexact_si(igusa_psi4(I), &I[0], 4);
  fmpz_divexact_si(igusa_psi6(I), &I[1], 4);
  fmpz_divexact_si(igusa_chi10(I), &I[2], -n_pow(2,12));
  fmpz_divexact_si(igusa_chi12(I), &I[3], n_pow(2,15));
  cov_normalize_fmpz(I, I, 4, weights);  
}

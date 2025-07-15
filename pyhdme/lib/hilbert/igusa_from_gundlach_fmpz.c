
#include "hilbert.h"

void igusa_from_gundlach_fmpz(fmpz* I, fmpz* G, slong delta)
{
  fmpz_t temp;
  slong weights[4] = IGUSA_HALFWEIGHTS;
  
  fmpz_init(temp);
  
  if (delta != 5)
    {
      flint_printf("(igusa_from_gundlach_fmpz) Error: Gundlach invariants only implemented for discriminant 5\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_mul(igusa_psi4(I), &G[0], &G[0]);

  fmpz_pow_ui(igusa_psi6(I), &G[0], 3);
  fmpz_mul_si(temp, &G[1], - 32*27);
  fmpz_add(igusa_psi6(I), igusa_psi6(I), temp);

  fmpz_neg(igusa_chi10(I), &G[2]);

  fmpz_mul(igusa_chi12(I), &G[1], &G[1]);
  fmpz_mul_si(igusa_chi12(I), igusa_chi12(I), 3);
  fmpz_mul(temp, &G[0], &G[2]);
  fmpz_mul_si(temp, temp, -2);
  fmpz_add(igusa_chi12(I), igusa_chi12(I), temp);

  cov_normalize_fmpz(I, I, 4, weights);
  
  fmpz_clear(temp);
}

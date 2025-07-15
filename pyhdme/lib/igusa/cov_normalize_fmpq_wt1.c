
#include "igusa.h"

void cov_normalize_fmpq_wt1(fmpz* I, fmpq* S, slong nb)
{
  fmpz_t den;
  fmpq_t num;
  slong k;

  fmpz_init(den);
  fmpq_init(num);

  fmpz_one(den);
  for (k = 0; k < nb; k++) fmpz_lcm(den, den, fmpq_denref(&S[k]));
  for (k = 0; k < nb; k++)
    {
      fmpq_mul_fmpz(num, &S[k], den);
      fmpq_numerator(&I[k], num);
    }
  cov_normalize_fmpz_wt1(I, I, nb);

  fmpz_clear(den);
  fmpq_clear(num);
}

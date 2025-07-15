
#include "igusa.h"

static int
igusa_I6prime_from_IC_fmpz(fmpz_t I6prime, fmpz* I)
{
  fmpz_t res, temp;
  int r;

  fmpz_init(res);
  fmpz_init(temp);
  
  /* Get I6prime from I6 */
  fmpz_mul_si(res, &I[2], -3);
  fmpz_mul(temp, &I[0], &I[1]);
  fmpz_add(res, res, temp);
  r = fmpz_divisible_si(res, 2);
  if (r)
    {
      fmpz_divexact_si(res, res, 2);
      fmpz_set(I6prime, res);
    }
  
  fmpz_clear(res);
  fmpz_clear(temp);
  return r;
}

void igusa_from_IC_fmpz(fmpz* I, fmpz* IC)
{
  fmpz_t I6prime, I12;
  slong weights2[4] = IC_HALFWEIGHTS;
  fmpz* resc;
  fmpz* S;

  fmpz_init(I6prime);
  fmpz_init(I12);
  resc = _fmpz_vec_init(4);
  S = _fmpz_vec_init(4);
  
  _fmpz_vec_set(resc, IC, 4);

  cov_rescale_fmpz_si(resc, IC, 2, 4, weights2);
  igusa_I6prime_from_IC_fmpz(I6prime, resc);
  fmpz_mul(I12, &resc[0], &resc[3]);

  fmpz_set(&S[0], &resc[1]);
  fmpz_set(&S[1], I6prime);
  fmpz_set(&S[2], &resc[3]);
  fmpz_set(&S[3], I12);

  igusa_from_streng_fmpz(I, S);

  fmpz_clear(I6prime);
  fmpz_clear(I12);
  _fmpz_vec_clear(resc, 4);
  _fmpz_vec_clear(S, 4);
}

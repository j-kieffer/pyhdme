
#include "igusa.h"

static int
igusa_I2_from_streng_fmpz(fmpz_t I2, fmpz* S)
{
  if (!fmpz_divisible(&S[3], &S[2]))
    {
      return 0;
    }
  else
    {
      fmpz_divexact(I2, &S[3], &S[2]);
      return 1;
    }      
}

static int
igusa_I6_from_streng_fmpz(fmpz_t I6, fmpz* S)
{
  fmpz_t I2, res, temp;
  int r;

  fmpz_init(I2);
  fmpz_init(res);
  fmpz_init(temp);

  r = igusa_I2_from_streng_fmpz(I2, S);
  if (!r) return 0;
  
  /* Get I6 from I6prime */
  fmpz_mul_si(res, &S[1], 2);
  fmpz_mul(temp, I2, &S[0]);
  fmpz_sub(res, res, temp);
  r = fmpz_divisible_si(res, 3);
  if (r)
    {
      fmpz_divexact_si(res, res, -3);
      fmpz_set(I6, res);
    }

  fmpz_clear(I2);
  fmpz_clear(res);
  fmpz_clear(temp);
  return r;
}


void igusa_IC_fmpz(fmpz* IC, fmpz* I)
{
  fmpz_t I2, I6;
  fmpz* S;
  fmpz* resc;
  slong weights[4] = IGUSA_WEIGHTS;
  slong weights2[4] = IC_HALFWEIGHTS;
  int r;

  fmpz_init(I2);
  fmpz_init(I6);
  resc = _fmpz_vec_init(4);
  S = _fmpz_vec_init(4);

  igusa_streng_fmpz(S, I);  
  cov_rescale_fmpz(resc, S, &S[2], 4, weights); /* Ensure I2 is an integer */
  cov_rescale_fmpz_si(resc, resc, 3, 4, weights); /* Ensure I6 is an integer */
  
  r = igusa_I2_from_streng_fmpz(I2, resc)
    && igusa_I6_from_streng_fmpz(I6, resc);
  if (!r)
    {
      flint_printf("(igusa_IC_fmpz) Error: I2 or I6 is not integral\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_set(&IC[0], I2);
  fmpz_set(&IC[1], &resc[0]);
  fmpz_set(&IC[2], I6);
  fmpz_set(&IC[3], &resc[2]);
  cov_normalize_fmpz(IC, IC, 4, weights2);

  fmpz_clear(I2);
  fmpz_clear(I6);
  _fmpz_vec_clear(resc, 4);
  _fmpz_vec_clear(S, 4);
}


#include "igusa.h"

void igusa_make_integral(fmpq_t scal, fmpz* I, slong wt, slong corr_2_3)
{
  fmpq* X;
  fmpz_factor_t fac;
  fmpz_t p;
  fmpz_t c;
  fmpq_t val;
  fmpq_t corr;
  slong k;
  slong weights[X_NB] = X_WEIGHTS;
  
  X = _fmpq_vec_init(X_NB);
  fmpz_factor_init(fac);
  fmpz_init(p);
  fmpz_init(c);
  fmpq_init(val);
  fmpq_init(corr);
  
  fmpz_one(p);
  igusa_X(X, I);
  cov_factors(fac, I, 4);
  
  fmpq_one(scal);
  for (k = 0; k < 2 + cov_factor_nb(fac); k++)
    {
      if (k == 0) fmpz_set_si(p, 2);
      else if (k == 1) fmpz_set_si(p, 3);
      else
	{
	  fmpz_set(p, cov_factor_p(fac, k-2));
	  if (fmpz_equal_si(p, 2) || fmpz_equal_si(p, 3)) continue;
	}
      
      cov_valuation_fmpq(val, X, p, X_NB, weights);

      /* Correction from superspecial reduction */
      if (fmpz_equal_si(p,2) || fmpz_equal_si(p,3))
	{
	  fmpq_mul_si(val, val, wt - corr_2_3);
	}
      else
	{
	  fmpq_mul_si(val, val, wt);
	}
      
      /* Get ceil(val) */
      fmpz_cdiv_q(c, fmpq_numref(val), fmpq_denref(val));
      if (fmpz_cmp_si(c, 0) < 0)
	{
	  fmpz_neg(c, c);
	  fmpz_pow_fmpz(c, p, c);
	  fmpq_mul_fmpz(scal, scal, c);
	}
      else
	{
	  fmpz_pow_fmpz(c, p, c);
	  fmpq_div_fmpz(scal, scal, c);
	}
    }

  _fmpq_vec_clear(X, X_NB);
  fmpz_factor_clear(fac);
  fmpz_clear(p);
  fmpz_clear(c);
  fmpq_clear(val);
  fmpq_clear(corr);
}

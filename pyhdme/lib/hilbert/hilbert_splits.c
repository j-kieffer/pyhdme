
#include "hilbert.h"

/* A very naive algorithm using exhaustive search up to some reasonable bound */
int hilbert_splits(fmpz_poly_t beta, slong ell, slong delta)
{
  fmpz_t a, b, discr;
  int stop = 0;
  int res;
  
  fmpz_init(a);
  fmpz_init(b);
  fmpz_init(discr);

  res = n_sqrtmod(delta % ell, ell);
  if (res == 0) stop = 1;
    
  fmpz_poly_zero(beta);
  while (!stop && fmpz_cmp_si(a, n_pow(10, 5)) < 0)
    {
      fmpz_add_si(a, a, 1);
      if (delta%2 == 0) /* Solve equation a^2 - b^2 delta/4 = ell */
	{
	  fmpz_mul(discr, a, a);
	  fmpz_sub_si(discr, discr, ell);
	  fmpz_mul_si(discr, discr, 4);
	  if (fmpz_fdiv_ui(discr, delta) == 0)
	    {
	      fmpz_divexact_si(discr, discr, delta);
	      if (fmpz_is_square(discr))
		{
		  stop = 1;
		  fmpz_sqrt(b, discr);
		  fmpz_poly_set_coeff_fmpz(beta, 0, a);
		  fmpz_poly_set_coeff_fmpz(beta, 1, b);
		}
	    }
	}
      else /* delta is odd, solve equation a^2 + ab + b^2(1-delta)/4 = ell */
	{
	  fmpz_mul(discr, a, a);
	  fmpz_sub_si(discr, discr, ell);
	  fmpz_mul_si(discr, discr, delta);
	  fmpz_add_si(discr, discr, ell);
	  fmpz_mul_si(discr, discr, 4);
	  if (fmpz_is_square(discr))
	    {
	      fmpz_sqrt(discr, discr);
	      fmpz_submul_ui(discr, a, 2);
	      if (fmpz_fdiv_ui(discr, delta-1) == 0)
		{
		  stop = 1;
		  fmpz_divexact_si(b, discr, 1-delta);
		  fmpz_poly_set_coeff_fmpz(beta, 0, a);
		  fmpz_poly_set_coeff_fmpz(beta, 1, b);
		}
	    }
	}
    } 
  
  res = hilbert_is_totally_positive(beta, delta);

  fmpz_clear(a);
  fmpz_clear(b);
  fmpz_clear(discr);
  return res;
}

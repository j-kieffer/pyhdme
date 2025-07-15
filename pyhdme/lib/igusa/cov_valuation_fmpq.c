
#include "igusa.h"

void cov_valuation_fmpq(fmpq_t val, fmpq* I, const fmpz_t p, slong nb, slong* weights)
{
  fmpq_t v;
  fmpz_t r;
  slong n1, n2;
  slong k;  
  int set = 0;

  fmpq_init(v);
  fmpz_init(r);
  
  for (k = 0; k < nb; k++)
    {
      if (!fmpq_is_zero(&I[k]))
	{
	  n1 = fmpz_remove(r, fmpq_numref(&I[k]), p);
	  n2 = fmpz_remove(r, fmpq_denref(&I[k]), p);
	  fmpq_set_si(v, n1 - n2, weights[k]);
	  if (!set || fmpq_cmp(v, val) < 0)
	    {
	      fmpq_set(val, v);
	      set = 1;
	    }
	}
    }
  if (!set)
    {
      flint_printf("(cov_valuation_fmpq) Error: all entries are zero\n");
      fflush(stdout);
      flint_abort();
    }

  fmpq_clear(v);
  fmpz_clear(r);  
}

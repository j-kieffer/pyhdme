
#include "hecke.h"

void siegel_T1_coset(fmpz_mat_t m, slong k, slong p)
{
  slong a, b, c;
  slong i;
  
  if ((k < 0) || (k >= siegel_nb_T1_cosets(p)))
    {
      flint_printf("(siegel_T1_coset) Error: no matrix numbered %wd\n", k);
      fflush(stdout);
      flint_abort();
    }

  fmpz_mat_zero(m);
  
  if (k == 0)
    {
      /* Case 1 */
      fmpz_set_si(fmpz_mat_entry(m, 0, 0), p);
      fmpz_set_si(fmpz_mat_entry(m, 1, 1), n_pow(p, 2));
      fmpz_set_si(fmpz_mat_entry(m, 2, 2), p);
      fmpz_set_si(fmpz_mat_entry(m, 3, 3), 1);
    }
  else if (k < 1 + (n_pow(p, 2)-1) )
    {
      /* Case 2 */
      if (k < 1 + (p-1))
	{
	  /* a is zero, b too, c is anything nonzero */
	  a = 0;
	  b = 0;
	  c = k;
	}
      else
	{
	  /* a is nonzero, b is anything, c is b^2/a */
	  /* k-p is between 0 and p(p-1)-1 */
	  a = (k-p) % (p-1);
	  a += 1;
	  b = (k-p) % p;
	  c = (b*b) % p;
	  c *= n_invmod(a, p);
	  c = c % p;
	}
      for (i = 0; i < 4; i++) fmpz_set_si(fmpz_mat_entry(m, i, i), p);
      fmpz_set_si(fmpz_mat_entry(m, 0, 2), a);
      fmpz_set_si(fmpz_mat_entry(m, 0, 3), b);
      fmpz_set_si(fmpz_mat_entry(m, 1, 2), b);
      fmpz_set_si(fmpz_mat_entry(m, 1, 3), c);
    }
  else if (k < n_pow(p, 2) + p)
    {
      /* Case 3 */
      a = k - n_pow(p, 2);
      fmpz_set_si(fmpz_mat_entry(m, 0, 0), n_pow(p, 2));
      fmpz_set_si(fmpz_mat_entry(m, 1, 0), -a*p);
      fmpz_set_si(fmpz_mat_entry(m, 1, 1), p);
      fmpz_set_si(fmpz_mat_entry(m, 2, 2), 1);
      fmpz_set_si(fmpz_mat_entry(m, 2, 3), a);
      fmpz_set_si(fmpz_mat_entry(m, 3, 3), p);
    }
  else if (k < n_pow(p, 2) + p + n_pow(p, 3))
    {
      /* Case 4 */
      k = k - n_pow(p, 2) - p;      
      b = k % p;
      a = k / p;
      fmpz_set_si(fmpz_mat_entry(m, 0, 0), 1);
      fmpz_set_si(fmpz_mat_entry(m, 0, 2), a);
      fmpz_set_si(fmpz_mat_entry(m, 0, 3), -b);
      fmpz_set_si(fmpz_mat_entry(m, 1, 1), p);
      fmpz_set_si(fmpz_mat_entry(m, 1, 2), -p*b);
      fmpz_set_si(fmpz_mat_entry(m, 2, 2), n_pow(p, 2));
      fmpz_set_si(fmpz_mat_entry(m, 3, 3), p);
    }
  else
    {
      /* Case 5 */
      k = k - n_pow(p, 3) - n_pow(p, 2) - p;
      a = k%p;
      k = k/p;
      b = k%p;
      c = k/p;
      fmpz_set_si(fmpz_mat_entry(m, 0, 0), p);
      fmpz_set_si(fmpz_mat_entry(m, 0, 3), b*p);
      fmpz_set_si(fmpz_mat_entry(m, 1, 0), -a);
      fmpz_set_si(fmpz_mat_entry(m, 1, 1), 1);
      fmpz_set_si(fmpz_mat_entry(m, 1, 2), b);
      fmpz_set_si(fmpz_mat_entry(m, 1, 3), a*b+c);
      fmpz_set_si(fmpz_mat_entry(m, 2, 2), p);
      fmpz_set_si(fmpz_mat_entry(m, 2, 3), a*p);
      fmpz_set_si(fmpz_mat_entry(m, 3, 3), n_pow(p, 2));	 
    }
      
}


#include "hecke.h"

void siegel_coset(fmpz_mat_t m, slong k, slong ell)
{
  slong a = 0;
  slong b = 0;
  slong c = 0;
  fmpz_mat_t eta, etaR;
  fmpz_mat_t temp;
  slong u, v, w, x, y, z;
  
  fmpz_mat_init(eta, 4, 4);
  fmpz_mat_init(etaR, 4, 4);
  fmpz_mat_init(temp, 4, 4);
  
  if ((k < 0) || (k >= siegel_nb_cosets(ell)))
    {
      flint_printf("(siegel_coset) No matrix numbered %wd\n", k);
      fflush(stdout);
      flint_abort();
    }
  
  if (k < n_pow(ell, 3))
    {
      /* Case T_1(a, b, c) */
      a = k % ell;
      k = k / ell;
      b = k % ell;
      k = k / ell;
      c = k % ell;
      
      fmpz_mat_one(temp);
      fmpz_mat_neg(temp, temp);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 2), a);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 3), b);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 2), b);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 3), c);
      fmpz_mat_set(eta, temp);
      
      fmpz_mat_one(etaR);
    }
  else if (k < n_pow(ell, 3) + ell + (ell - 1) * ell)
    {
      /* Case T_2(a, b, c) */
      if (k < n_pow(ell, 3) + ell)
	{
	  a = 0;
	  b = 0;
	  c = k - n_pow(ell, 3);

	  /* Special etaR in this case */
	  if (c == 0)
	    {
	      fmpz_mat_J(etaR);
	    }
	  else
	    {
	      n_xgcd((ulong*) &v, (ulong*) &u, ell, c); 
	      u = -u; /* uc + lv = 1 */
	      fmpz_mat_zero(temp);
	      fmpz_set_si(fmpz_mat_entry(temp, 0, 1), u);
	      fmpz_set_si(fmpz_mat_entry(temp, 0, 3), v);
	      fmpz_set_si(fmpz_mat_entry(temp, 1, 2), 1);
	      fmpz_set_si(fmpz_mat_entry(temp, 2, 1), -ell);
	      fmpz_set_si(fmpz_mat_entry(temp, 2, 3), c);
	      fmpz_set_si(fmpz_mat_entry(temp, 3, 0), -1);
	      fmpz_mat_set(etaR, temp);
	    }
	}
      else
	{
	  a = (k - n_pow(ell, 3) - ell) % (ell-1);
	  a += 1;
	  b = (k - n_pow(ell, 3) - ell) % ell;
	  c = (b*b) % ell;
	  c *= n_invmod(a, ell);
	  c = c % ell;

	  n_xgcd((ulong*) &v, (ulong*) &u, ell, a);
	  u = -u; /* au + lv = 1 */
	  w = - b*u;
	  w = w % ell; /* w = -b/a mod ell, u = 1/a mod ell */
	  x = - (b + w*a) / ell;
	  y = - (c + w*b) / ell;
	  z = u*x - w*v;

	  fmpz_mat_zero(temp);
	  fmpz_set_si(fmpz_mat_entry(temp, 0, 0), u);
	  fmpz_set_si(fmpz_mat_entry(temp, 0, 2), v);
	  fmpz_set_si(fmpz_mat_entry(temp, 0, 3), z);
	  fmpz_set_si(fmpz_mat_entry(temp, 1, 3), -1);
	  fmpz_set_si(fmpz_mat_entry(temp, 2, 0), -ell);
	  fmpz_set_si(fmpz_mat_entry(temp, 2, 2), a);
	  fmpz_set_si(fmpz_mat_entry(temp, 2, 3), b);
	  fmpz_set_si(fmpz_mat_entry(temp, 3, 0), w);
	  fmpz_one(fmpz_mat_entry(temp, 3, 1));
	  fmpz_set_si(fmpz_mat_entry(temp, 3, 2), x);
	  fmpz_set_si(fmpz_mat_entry(temp, 3, 3), y);
	  fmpz_mat_set(etaR, temp);
	}
      /* eta is always given by the same formula */
      fmpz_mat_zero(temp);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 0), -a);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 1), -b);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 0), -b);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 1), -c);
      fmpz_one(fmpz_mat_entry(temp, 0, 2));
      fmpz_one(fmpz_mat_entry(temp, 1, 3));
      fmpz_set_si(fmpz_mat_entry(temp, 2, 0), -1);
      fmpz_set_si(fmpz_mat_entry(temp, 3, 1), -1);
      fmpz_mat_set(eta, temp);
    }
  else if (k < n_pow(ell, 3) + ell*ell + ell)
    {
      /* Case T_3(a) */
      a = k - n_pow(ell, 3) - ell*ell;

      fmpz_mat_zero(temp);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 0), -1);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 1), -a);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 2), -a);
      fmpz_one(fmpz_mat_entry(temp, 1, 3));
      fmpz_set_si(fmpz_mat_entry(temp, 2, 2), -1);
      fmpz_set_si(fmpz_mat_entry(temp, 3, 1), -1);
      fmpz_mat_set(eta, temp);

      fmpz_mat_zero(temp);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 3), -1);
      fmpz_one(fmpz_mat_entry(temp, 1, 0));
      fmpz_one(fmpz_mat_entry(temp, 2, 1));
      fmpz_one(fmpz_mat_entry(temp, 3, 2));
      fmpz_mat_set(etaR, temp);
    }
  else
    {
      /* Case T_4 */
      fmpz_mat_zero(temp);
      fmpz_one(fmpz_mat_entry(temp, 0, 1));
      fmpz_one(fmpz_mat_entry(temp, 1, 1));
      fmpz_one(fmpz_mat_entry(temp, 1, 2));
      fmpz_one(fmpz_mat_entry(temp, 2, 0));
      fmpz_set_si(fmpz_mat_entry(temp, 2, 1), -1);
      fmpz_one(fmpz_mat_entry(temp, 2, 2));
      fmpz_one(fmpz_mat_entry(temp, 2, 3));
      fmpz_set_si(fmpz_mat_entry(temp, 3, 0), -1);
      fmpz_one(fmpz_mat_entry(temp, 3, 1));
      fmpz_mat_set(eta, temp);

      fmpz_mat_zero(temp);
      fmpz_one(fmpz_mat_entry(temp, 0, 0));
      fmpz_set_si(fmpz_mat_entry(temp, 1, 1), ell);
      fmpz_one(fmpz_mat_entry(temp, 1, 3));
      fmpz_set_si(fmpz_mat_entry(temp, 2, 0), -ell);
      fmpz_set_si(fmpz_mat_entry(temp, 2, 1), ell);
      fmpz_one(fmpz_mat_entry(temp, 2, 2));
      fmpz_one(fmpz_mat_entry(temp, 2, 3));
      fmpz_one(fmpz_mat_entry(temp, 3, 0));
      fmpz_set_si(fmpz_mat_entry(temp, 3, 1), -1);
      fmpz_mat_set(etaR, temp);
    }
  
  if (!fmpz_mat_is_symplectic(eta) || !fmpz_mat_is_symplectic(etaR))
    {
      fmpz_mat_print(eta);
      fmpz_mat_print(etaR);
      flint_printf("k = %wd, a = %wd, b = %wd, c = %wd, ell = %wd\n", k, a, b, c, ell);
      flint_printf("eta is symplectic? %wd; etaR is symplectic? %wd\n",
		   fmpz_mat_is_symplectic(eta), fmpz_mat_is_symplectic(etaR));
      fflush(stdout);
      flint_abort();
    }
  /* Compute transformation directly: warning, they are only general
     symplectic */
  for (u = 2; u < 4; u++)
    {
      for (v = 0; v < 4; v++)
	{
	  fmpz_mul_si(fmpz_mat_entry(eta, u, v),
		      fmpz_mat_entry(eta, u, v), ell);
	}
    }
  fmpz_mat_mul(m, etaR, eta);
  
  fmpz_mat_clear(eta);
  fmpz_mat_clear(etaR);
  fmpz_mat_clear(temp);
}

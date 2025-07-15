
#include "siegel.h"

static int
arb_round_fmpz(fmpz_t z, const arb_t x, const arb_t tol)
{
  arb_t rad;
  int res;
  arb_init(rad);
  
  res = 1;
  arb_get_rad_arb(rad, x);
  
  if (!arb_lt(rad, tol)) res = 0;
  if (res)
    {
      arf_get_fmpz(z, arb_midref(x), ARF_RND_NEAR);
    }
  
  arb_clear(rad);
  return res;
}


int siegel_reduce_real(acb_mat_t w, fmpz_mat_t u, const acb_mat_t z,
		       const arb_t tol, slong prec)
{
  int res = 1;
  fmpz_mat_t round, zero, one;
  slong g = arb_mat_nrows(z);
  slong i, j;

  fmpz_mat_init(round, g, g);
  fmpz_mat_init(zero, g, g);
  fmpz_mat_init(one, g, g);
  fmpz_mat_one(one);
  
  for (i = 0; i < g; i++)
    {
      for (j = 0; j < i; j++)
	{
	  fmpz_set(acb_mat_entry(round, i, j),
		   acb_mat_entry(round, j, i));
	}
      for (j = i; j < g; j++)
	{
	  res = res && arb_round_fmpz(fmpz_mat_entry(round, i, j),
				      acb_realref(acb_mat_entry(z, i, j)), tol);
	}
    }
  fmpz_mat_neg(round, round);
  fmpz_mat_set_abcd(u, one, round, zero, one);
  siegel_transform(w, u, z, prec);

  fmpz_mat_clear(round);
  fmpz_mat_clear(zero);
  fmpz_mat_clear(one);
  return res;
}

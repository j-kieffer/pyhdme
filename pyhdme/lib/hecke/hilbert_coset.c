
#include "hecke.h"

void hilbert_coset(fmpz_poly_mat_t m, slong k,
		   const fmpz_poly_t beta, slong ell, slong delta)
{
  if (k == 0)
    {
      fmpz_poly_zero(fmpz_poly_mat_entry(m, 0, 0));
      fmpz_poly_one(fmpz_poly_mat_entry(m, 0, 1));
      fmpz_poly_set_si(fmpz_poly_mat_entry(m, 1, 0), -1);
      fmpz_poly_zero(fmpz_poly_mat_entry(m, 1, 1));
    }
  else
    {
      fmpz_poly_one(fmpz_poly_mat_entry(m, 0, 0));
      fmpz_poly_set_si(fmpz_poly_mat_entry(m, 0, 1), k % ell);
      fmpz_poly_zero(fmpz_poly_mat_entry(m, 1, 0));
      fmpz_poly_one(fmpz_poly_mat_entry(m, 1, 1));
    }
  
  fmpz_poly_mul(fmpz_poly_mat_entry(m, 1, 0),
		fmpz_poly_mat_entry(m, 1, 0), beta);
  fmpz_poly_mul(fmpz_poly_mat_entry(m, 1, 1),
		fmpz_poly_mat_entry(m, 1, 1), beta);
}

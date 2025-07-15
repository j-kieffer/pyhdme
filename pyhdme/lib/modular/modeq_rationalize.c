
#include "modular.h"

int modeq_rationalize(modeq_t R, const modeq_acb_t E, slong prec)
{
  fmpz* den_vec;
  acb_poly_struct* quo_vec;
  fmpq_poly_struct* pol_vec;
  slong k;
  slong nb = 1+modeq_nb(E);
  int res = 1;

  den_vec = flint_malloc(nb * sizeof(fmpz));
  pol_vec = flint_malloc(nb * sizeof(fmpq_poly_struct));
  quo_vec = flint_malloc(nb * sizeof(acb_poly_struct));
  for (k = 0; k < nb; k++)
    {
      fmpz_init(&den_vec[k]);
      fmpq_poly_init(&pol_vec[k]);
      acb_poly_init(&quo_vec[k]);
    }

  /* Rationalize input */
  for (k = 0; k < nb; k++)
    {
      if (res)
	{
	  acb_poly_scalar_div(&quo_vec[k], &modeq_all_nums(E)[k], modeq_den(E), prec);
	  res = acb_poly_rationalize(&pol_vec[k], &quo_vec[k],
				     modeq_degree(E), prec);
	}
    }

  /* Set R */
  modeq_degree(R) = modeq_degree(E);
  modeq_nb(R) = modeq_nb(E);

  if (res)
    { 
      for (k = 0; k < nb; k++)
	{
	  fmpq_poly_get_denominator(&den_vec[k], &pol_vec[k]);
	}
      _fmpz_vec_lcm(modeq_den(R), den_vec, nb);      
      for (k = 0; k < nb; k++)
	{
	  fmpq_poly_scalar_mul_fmpz(&pol_vec[k], &pol_vec[k], modeq_den(R));
	  fmpq_poly_get_numerator(&modeq_all_nums(R)[k], &pol_vec[k]);
	}
    }
  if (res) pol_simplify(modeq_all_nums(R), modeq_den(R), modeq_degree(R),
			modeq_nb(R)+1);

  for (k = 0; k < nb; k++)
    {
      fmpz_clear(&den_vec[k]);
      fmpq_poly_clear(&pol_vec[k]);
      acb_poly_clear(&quo_vec[k]);
    }
  flint_free(den_vec);
  flint_free(pol_vec);
  flint_free(quo_vec);
  return res;
}

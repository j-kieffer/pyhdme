
#include "modular.h"

void modeq_isog_monomials_nf(fmpz_poly_struct* M, const modeq_t E,
			     const fmpz_poly_t factor, slong mult)
{  
  fmpz_poly_t num;
  fmpq* aux;
  slong nb = modeq_nb(E);
  slong j;
  
  fmpz_poly_init(num);
  aux = _fmpq_vec_init(nb);

  /* Evaluate as number field elements */
  for (j = 0; j < modeq_nb(E); j++)
    {
      pol_remove_factor_Q(num, modeq_interpolate(E, j), factor, mult-1);
      fmpz_poly_rem(&M[j], num, factor); /* Might be higher degree, but still integer coefs? */
    }
  
  fmpz_poly_clear(num);
  _fmpq_vec_clear(aux, nb);
}


#include "modular.h"

void modeq_isog_monomials_Q(fmpz* M, const modeq_t E, const fmpq_t root, slong mult)
{
  fmpz_poly_t num;
  fmpq* aux;
  slong nb = modeq_nb(E);
  slong j;
  
  fmpz_poly_init(num);
  aux = _fmpq_vec_init(nb);

  /* Evaluate as rational numbers */
  for (j = 0; j < modeq_nb(E); j++)
    {
      fmpz_poly_set(num, modeq_interpolate(E, j));
      pol_remove_root_Q(num, num, root, mult-1);
      fmpz_poly_evaluate_fmpq(&aux[j], num, root);
    }

  /* Rescale to integers */
  cov_normalize_fmpq_wt1(M, aux, nb);

  fmpz_poly_clear(num);
  _fmpq_vec_clear(aux, nb);
}


#include "modular.h"

void modeq_clear(modeq_t E)
{
  slong k;
  slong nb = MODEQ_MAX_NB_MONOMIALS;
  
  fmpz_clear(modeq_den(E));
  fmpz_poly_clear(modeq_equation(E));
  for (k = 0; k < nb; k++) fmpz_poly_clear(modeq_interpolate(E, k));
  flint_free(modeq_all_nums(E));
}

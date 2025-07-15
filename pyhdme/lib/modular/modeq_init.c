
#include "modular.h"

void modeq_init(modeq_t E)
{
  slong k;
  slong nb = MODEQ_MAX_NB_MONOMIALS;
  
  fmpz_init(modeq_den(E));  
  modeq_all_nums(E) = flint_malloc((nb+1) * sizeof(fmpz_poly_struct));
  fmpz_poly_init(modeq_equation(E));
  for (k = 0; k < nb; k++) fmpz_poly_init(modeq_interpolate(E, k));
}

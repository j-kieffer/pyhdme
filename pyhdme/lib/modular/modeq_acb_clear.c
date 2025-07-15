
#include "modular.h"

void modeq_acb_clear(modeq_acb_t E)
{
  slong k;
  slong nb = MODEQ_MAX_NB_MONOMIALS;
  
  acb_clear(modeq_den(E));
  acb_poly_clear(modeq_equation(E));
  for (k = 0; k < nb; k++) acb_poly_clear(modeq_interpolate(E, k));
  flint_free(modeq_all_nums(E));
}

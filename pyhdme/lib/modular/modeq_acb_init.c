
#include "modular.h"

void modeq_acb_init(modeq_acb_t E)
{
  slong k;
  slong nb = MODEQ_MAX_NB_MONOMIALS;
  
  acb_init(modeq_den(E));
  modeq_all_nums(E) = flint_malloc((nb+1) * sizeof(acb_poly_struct));
  acb_poly_init(modeq_equation(E));
  for (k = 0; k < nb; k++) acb_poly_init(modeq_interpolate(E, k));
}

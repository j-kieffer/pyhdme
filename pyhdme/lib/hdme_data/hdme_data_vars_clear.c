
#include "hdme_data.h"

void hdme_data_vars_clear(char** vars, slong nb)
{
  slong k;
  for (k = 0; k < nb; k++) flint_free(vars[k]);
  flint_free(vars);
}


#include "hdme_data.h"

char** hdme_data_vars_init(slong nb)
{
  char** vars;
  slong k;

  vars = flint_malloc(nb * sizeof(char*));
  for (k = 0; k < nb; k++)
    {
      vars[k] = flint_malloc((HDME_DATA_VAR_LEN + 1) * sizeof(char));
    }
  return vars;
}

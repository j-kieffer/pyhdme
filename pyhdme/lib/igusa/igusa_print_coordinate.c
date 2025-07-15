
#include "igusa.h"

void igusa_print_coordinate(const fmpz_mpoly_t pol, const fmpz_mpoly_ctx_t ctx)
{
  char** x;
  x = hdme_data_vars_init(4);
  hdme_data_vars_set(x, "psi4", 0);
  hdme_data_vars_set(x, "psi6", 1);
  hdme_data_vars_set(x, "chi10", 2);
  hdme_data_vars_set(x, "chi12", 3);  
  fmpz_mpoly_print_pretty(pol, (const char**) x, ctx);
  hdme_data_vars_clear(x, 4);
}

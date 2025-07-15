
#include <stdint.h>
#include <string.h>

#include "hdme_data.h"

void hdme_data_read(fmpq_mpoly_t pol, const char** vars, const char* name,
		    const fmpq_mpoly_ctx_t ctx)
{
  char str[HDME_DATA_STR_LEN];
  int res;

  hdme_data_get_string(str, name);
  res = fmpq_mpoly_set_str_pretty(pol, str, vars, ctx);
  if (res == -1)
    {
      flint_printf("(hdme_data_read) Error parsing string: %s\n", str);
      flint_printf("(hdme_data_read) Key name: %s\n", name);
      fflush(stdout);
      flint_abort();
    }
}

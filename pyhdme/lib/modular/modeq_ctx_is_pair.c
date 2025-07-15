
#include "modular.h"

int modeq_ctx_is_pair(slong i1, slong i2, const modeq_ctx_t ctx)
{
  slong k;
  int res = 0;
  
  for (k = 0; k < modeq_ctx_nb_pairs(ctx); k++)
    {
      if (modeq_ctx_pair(ctx, k)[0] == i1
	  && modeq_ctx_pair(ctx, k)[1] == i2)
	{
	  res = 1;
	  break;
	}
    }
  
  return res;
}

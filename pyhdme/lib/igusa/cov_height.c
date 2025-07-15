
#include "igusa.h"

slong cov_height(fmpz* I, slong len, slong* weights)
{
  slong k;
  slong h = 1;
  
  for (k = 0; k < len; k++)
    {
      h = FLINT_MAX(h, fmpz_bits(&I[k]) / weights[k]);
    }
  return h;
}

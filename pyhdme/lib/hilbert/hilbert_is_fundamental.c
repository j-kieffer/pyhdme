
#include "hilbert.h"

int hilbert_is_fundamental(slong delta)
{
  int res = 0;
  if ((delta >= 5) && (delta % 4 == 1) && n_is_squarefree(delta))
    {
      res = 1;
    }
  if ((delta >= 5) && (delta % 4 == 0) && n_is_squarefree(delta/4) && ((delta/4) % 4 != 1))
    {
      res = 1;
    }
  return res;
}

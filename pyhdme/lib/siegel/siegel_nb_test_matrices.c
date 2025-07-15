
#include "siegel.h"


slong siegel_nb_test_matrices(slong g)
{
  slong res;
  switch(g)
    {
    case 1:
      res = 1;
      break;
    case 2:
      res = 19;
      break;
    default:
      flint_printf("siegel_nb_test_matrices not implemented for g = %d\n", g);
      flint_abort();
    }
  return res;
}

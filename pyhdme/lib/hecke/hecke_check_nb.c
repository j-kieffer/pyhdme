
#include "hecke.h"


void hecke_check_nb(const hecke_t H, slong nb)
{
  if (hecke_nb(H) != nb)
    {
      flint_printf("(hecke_check_nb) Error: Hecke data structure initialized with %wd slots instead of the expected %wd\n", hecke_nb(H), nb);
      fflush(stdout);
      flint_abort();
    }
}

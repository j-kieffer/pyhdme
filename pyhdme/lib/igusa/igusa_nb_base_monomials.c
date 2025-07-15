
#include "igusa.h"

slong igusa_nb_base_monomials(slong w)
{
  if (w == 20) return 5;
  else if (w == 30) return 7;
  else if (w == 60) return 10;
  else
    {
      flint_printf("(cov_nb_base_monomials) Error: invalid weight %wd\n", w);
      fflush(stdout);
      flint_abort();
    }
  return 0;
}

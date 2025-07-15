
#include "igusa.h"

void igusa_base_exps(slong* exps, slong wt, slong k)
{
  slong j;
  slong* ptr;
  slong wt20[4*5] =
    {
      0, 0, 2, 0,
      1, 1, 1, 0,
      2, 0, 0, 1,
      5, 0, 0, 0,
      2, 2, 0, 0
    };  
  slong wt30[4*7] =
    {
      0, 5, 0, 0,
      3, 3, 0, 0,
      0, 0, 3, 0,
      5, 0, 1, 0,
      0, 3, 0, 1,
      1, 1, 2, 0,
      2, 0, 1, 1
    };
  slong wt60[4*10] =
    {
      0, 0, 6, 0,
      0, 0, 0, 5,
      0, 5, 3, 0,
      0, 2, 0, 4,
      3, 0, 0, 4,
      0, 10, 0, 0,
      3, 8, 0, 0,
      5, 0, 4, 0,
      1, 1, 5, 0,
      3, 2, 0, 3
    };

  if (k < 0 || k >= igusa_nb_base_monomials(wt))
    {
      flint_printf("(igusa_base_exps) Invalid index %wd\n", k);
      fflush(stdout);
      flint_abort();
    }
  if (wt == 20) ptr = &wt20[4*k];
  else if (wt == 30) ptr = &wt30[4*k];
  else if (wt == 60) ptr = &wt60[4*k];
  else
    {
      flint_printf("(igusa_base_exps) Invalid weight %wd\n", wt);
      fflush(stdout);
      flint_abort();
    }

  for (j = 0; j < 4; j++) exps[j] = ptr[j];    
}

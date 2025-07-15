
#include "igusa.h"

slong cov_nb_monomials(slong wt, slong nb, slong* weights)
{
  slong s;
  slong w;

  if (wt < 0)
    {
      flint_printf("(cov_nb_monomials) Error: negative target weight %wd\n", wt);
      fflush(stdout);
      flint_abort();
    }
  if (nb >= 1 && weights[0] <= 0)
    {
      flint_printf("(cov_nb_monomials) Error: nonpositive coordinate weight %wd\n",
		   weights[0]);
      fflush(stdout);
      flint_abort();
    }
  
  /* Base case: nb = 0 */
  if (nb == 0 && wt == 0) return 1;
  else if (nb == 0) return 0;

  /* Else: find out exponent of first variable and sum */
  w = 0;
  s = 0;
  while (w <= wt)
    {
      s += cov_nb_monomials(wt - w, nb - 1, &weights[1]);
      w += weights[0];
    }

  return s;
}

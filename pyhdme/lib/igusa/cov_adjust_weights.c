
#include "igusa.h"

void cov_adjust_weights(slong* adj, slong* weights, fmpz* I, slong nb)
{
  slong g = 0;
  slong k;
  
  for (k = 0; k < nb; k++)
    {
      if (!fmpz_is_zero(&I[k])) g = n_gcd(g, weights[k]);
    }
  
  if (g == 0)
    {
      flint_printf("(cov_adjust_weights) Error: all coordinates are zero\n");
      fflush(stdout);
      flint_abort();
    }
  
  for (k = 0; k < nb; k++)
    {
      /* We don't care about weight when I[k] = 0 */
      adj[k] = weights[k] / g;
    }
}

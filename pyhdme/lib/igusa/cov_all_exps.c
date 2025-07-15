
#include "igusa.h"

void cov_all_exps(slong* exps, slong wt, slong nb, slong* weights)
{
  slong k, i;
  slong cur_j, aux_j;
  slong* aux;
  slong w;

  if (nb >= 1 && weights[0] <= 0)
    {
      flint_printf("(cov_all_exps) Error: nonpositive coordinate weight\n");
      fflush(stdout);
      flint_abort();
    }
    
  if (wt == 0)
    {
      for (k = 0; k < nb; k++) exps[k] = 0;
    }
  
  else if (nb >= 1)
    {
      w = 0;
      cur_j = 0;
      while (w <= wt)
	{
	  aux_j = cov_nb_monomials(wt - w, nb - 1, &weights[1]);
	  aux = flint_malloc(aux_j * (nb-1) * sizeof(slong));
	  
	  cov_all_exps(aux, wt - w, nb - 1, &weights[1]);
	  for (k = cur_j; k < cur_j + aux_j; k++)
	    {
	      exps[k*nb] = w / weights[0];
	      for (i = 1; i < nb; i++)
		{
		  exps[k*nb + i] = aux[(k-cur_j)*(nb-1) + (i-1)];
		}
	    }

	  flint_free(aux);
	  cur_j += aux_j;
	  w += weights[0];	  
	}
    }  
}

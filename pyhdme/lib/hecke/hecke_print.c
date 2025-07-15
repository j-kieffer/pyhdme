
#include "hecke.h"

void hecke_print(const hecke_t H, slong digits)
{
  slong j, k;
  
  flint_printf("Hecke correspondence of level %wd\n", hecke_ell(H));
  flint_printf("Base matrix:\n");
  acb_mat_printd(hecke_tau(H), digits);
  flint_printf("Invariants:\n");
  for (k = 0; k < 4; k++)
    {
      acb_printd(&hecke_I_tau(H)[k], digits); flint_printf("\n");
    }
  flint_printf("%wd isogenous period matrices with cocycles and invariants:\n",
	       hecke_nb(H));
  for (j = 0; j < hecke_nb(H); j++)
    {
      acb_mat_printd(hecke_isog(H, j), digits);
      acb_mat_printd(hecke_star(H, j), digits);
      for (k = 0; k < 4; k++)
	{
	  acb_printd(&hecke_I(H, j)[k], digits); flint_printf("\n");
	}
    }  
}

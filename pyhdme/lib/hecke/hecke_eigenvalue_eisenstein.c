
#include "hecke.h"

void hecke_eigenvalue_eisenstein(fmpz_t eig, slong k, slong ell)
{
  fmpz_t ell_z;
  fmpz_t temp;

  fmpz_init(ell_z);
  fmpz_init(temp);
  
  if (k < 4)
    {
      flint_printf("(hecke_eigenvalue_eisenstein) No Eisenstein series of weight %wd\n", k);
      fflush(stdout);
      flint_abort();
    }

  fmpz_set_si(eig, 1);
  fmpz_set_si(ell_z, ell);
  
  fmpz_pow_ui(temp, ell_z, k-2);
  fmpz_add(eig, eig, temp);
  fmpz_pow_ui(temp, ell_z, k-1);
  fmpz_add(eig, eig, temp);
  fmpz_pow_ui(temp, ell_z, 2*k-3);
  fmpz_add(eig, eig, temp);

  fmpz_clear(ell_z);
  fmpz_clear(temp);
}


#include "hecke.h"

void hecke_eigenvalues_eisenstein_2(fmpz* eig, slong k, slong ell)
{
  fmpz_t ell_z;
  fmpz_t temp;
  fmpz_t lambda;

  fmpz_init(ell_z);
  fmpz_init(temp);
  fmpz_init(lambda);

  if (k < 4)
    {
      flint_printf("(hecke_eigenvalues_eisenstein_2) No Eisenstein series of weight %wd\n", k);
      fflush(stdout);
      flint_abort();
    }
  
  fmpz_set_si(ell_z, ell);

  /* lambda0 */
  fmpz_set_si(lambda, 1);
  
  fmpz_pow_ui(temp, ell_z, 4*k-6);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, 2*k-3);
  fmpz_addmul_ui(lambda, temp, 2);
  fmpz_pow_ui(temp, ell_z, 3*k-4);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, k-1);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, 2*k-2);
  fmpz_add(lambda, lambda, temp);
  
  fmpz_pow_ui(temp, ell_z, 2*k-4);
  fmpz_sub(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, 3*k-6);
  fmpz_sub(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, k-3);
  fmpz_sub(lambda, lambda, temp);
  
  fmpz_set(&eig[0], lambda);

  /* lambda1 */
  fmpz_zero(lambda);

  fmpz_pow_ui(temp, ell_z, 2*k-4);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, 2*k-6);
  fmpz_sub(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, 3*k-6);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, 3*k-5);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, k-3);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, k-2);
  fmpz_add(lambda, lambda, temp);

  fmpz_set(&eig[1], lambda);
  
  /* lambda2 */
  fmpz_pow_ui(&eig[2], ell_z, 2*k-6);

  /* lambda */
  fmpz_set_si(lambda, 1);

  fmpz_pow_ui(temp, ell_z, 4*k-6);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, k-1);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, k-1);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, 2*k-3);
  fmpz_addmul_ui(lambda, temp, 2);
  fmpz_pow_ui(temp, ell_z, 3*k-5);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, 3*k-4);
  fmpz_add(lambda, lambda, temp);
  fmpz_pow_ui(temp, ell_z, 2*k-2);
  fmpz_add(lambda, lambda, temp);

  fmpz_set(&eig[3], lambda);  

  fmpz_clear(ell_z);
  fmpz_clear(temp);
}

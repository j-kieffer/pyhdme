
#include "siegel.h"

void fmpz_mat_randtest_symplectic(fmpz_mat_t m, flint_rand_t state, slong bits)
{
  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t n;

  fmpz_mat_init(n, 2*g, 2*g);
  
  fmpz_mat_randtest_triangular_symplectic(m, state, bits);
  fmpz_mat_randtest_triangular_symplectic(m, state, bits);
  fmpz_mat_randtest_diagonal_symplectic(n, state, bits);
  fmpz_mat_mul(m, m, n);
  fmpz_mat_J(n);
  fmpz_mat_mul(m, m, n);
  fmpz_mat_randtest_triangular_symplectic(n, state, bits);
  fmpz_mat_mul(m, m, n);
  fmpz_mat_J(n);
  fmpz_mat_mul(m, m, n);
  fmpz_mat_randtest_diagonal_symplectic(n, state, bits);
  fmpz_mat_mul(m, m, n);
  fmpz_mat_J(n);
  fmpz_mat_mul(m, m, n);
  fmpz_mat_randtest_triangular_symplectic(n, state, bits);
  fmpz_mat_mul(m, m, n);

  fmpz_mat_clear(n);
}

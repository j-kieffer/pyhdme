
#include "siegel.h"

void fmpz_mat_randtest_diagonal_symplectic(fmpz_mat_t m, flint_rand_t state, slong bits)
{
  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t u;

  fmpz_mat_init(u, g, g);  
  bits = FLINT_MAX(bits, 1);
  
  fmpz_mat_one(u);
  fmpz_mat_randops(u, state, 2 * bits * g);
  fmpz_mat_diagonal_symplectic(m, u);

  fmpz_mat_clear(u);
}

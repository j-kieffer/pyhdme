
#include "siegel.h"

void fmpz_mat_randtest_triangular_symplectic(fmpz_mat_t m, flint_rand_t state, slong bits)
{
  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t zero, one;
  fmpz_mat_t b, bt;
  
  fmpz_mat_init(zero, g, g);
  fmpz_mat_init(one, g, g);
  fmpz_mat_init(b, g, g);
  fmpz_mat_init(bt, g, g);
  bits = FLINT_MAX(bits, 1);

  fmpz_mat_one(one);  
  fmpz_mat_randbits(b, state, bits);
  fmpz_mat_transpose(bt, b);
  fmpz_mat_add(b, b, bt);
  fmpz_mat_scalar_tdiv_q_2exp(b, b, 1);
  fmpz_mat_set_abcd(m, one, b, zero, one);

  fmpz_mat_clear(zero);
  fmpz_mat_clear(one);
  fmpz_mat_clear(b);
  fmpz_mat_clear(bt);
}

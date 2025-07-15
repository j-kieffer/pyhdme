
#include "modular.h"

void modeq_all_isog_Q(slong* nb_roots, fmpz* all_I,
		      const modeq_t E, const modeq_ctx_t ctx)
{
  slong nb = modeq_degree(E);
  slong nb_M = modeq_ctx_nb(ctx);
  slong* mults;
  fmpq* roots;
  fmpz* M;
  slong k;

  mults = flint_malloc(nb * sizeof(slong));
  roots = _fmpq_vec_init(nb);
  M = _fmpz_vec_init(nb_M);
  
  pol_roots_Q(nb_roots, roots, mults, modeq_equation(E));
  for (k = 0; k < *nb_roots; k++)
    {
      modeq_isog_monomials_Q(M, E, &roots[k], mults[k]);
      igusa_from_monomials(&all_I[4*k], M, modeq_ctx_weight(ctx));
    }

  flint_free(mults);
  _fmpq_vec_clear(roots, nb);
  _fmpz_vec_clear(M, nb_M);
}
  


#include "modular.h"

int modeq_all_isog_Fp(slong* nb_roots, fmpz* all_I, const modeq_t E,
		      const modeq_ctx_t ctx, const fmpz_mod_ctx_t fpctx)
{
  slong nb = modeq_degree(E);
  slong nb_M = modeq_ctx_nb(ctx);
  slong* mults;
  fmpz* roots;
  fmpz* M;
  fmpz_mod_poly_t red;
  slong k;
  int res;

  mults = flint_malloc(nb * sizeof(slong));
  roots = _fmpz_vec_init(nb);
  M = _fmpz_vec_init(nb_M);
  fmpz_mod_poly_init(red, fpctx);
  
  res = pol_reduce_Fp(red, modeq_equation(E), modeq_den(E), fpctx);
  if (res)
    {
      pol_roots_Fp(nb_roots, roots, mults, red, fpctx);
      for (k = 0; k < *nb_roots; k++)
	{
	  modeq_isog_monomials_Fp(M, E, &roots[k], mults[k], fpctx);
	  igusa_from_monomials(&all_I[4*k], M, modeq_ctx_weight(ctx));	  
	}
    }

  for (k = 0; k < 4*(*nb_roots); k++)
    {
      fmpz_mod(&all_I[k], &all_I[k], fmpz_mod_ctx_modulus(fpctx));
    }
  
  flint_free(mults);
  _fmpz_vec_clear(roots, nb);
  _fmpz_vec_clear(M, nb_M);
  fmpz_mod_poly_clear(red, fpctx);
  return res;
}

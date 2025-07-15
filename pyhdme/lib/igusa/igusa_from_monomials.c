
#include "igusa.h"

void igusa_from_monomials(fmpz* I, fmpz* M, slong wt)
{
  slong nb = igusa_nb_base_monomials(wt);
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t pol;
  int z4, z6, z10, z12;
  slong* e4;
  slong* e6;
  slong* e10;
  slong* e12;
  slong weights[4] = IGUSA_HALFWEIGHTS;
  fmpz* res;

  fmpz_mpoly_ctx_init(ctx, nb, ORD_LEX);
  fmpz_mpoly_init(pol, ctx);
  e4 = flint_malloc(nb * sizeof(slong));
  e6 = flint_malloc(nb * sizeof(slong));
  e10 = flint_malloc(nb * sizeof(slong));
  e12 = flint_malloc(nb * sizeof(slong));
  res = _fmpz_vec_init(4);
  
  igusa_from_monomials_zeroes(&z4, &z6, &z10, &z12, M, wt);
  igusa_from_monomials_exps(e4, e6, e10, e12, z4, z6, z10, z12, wt);

  cov_monomial(pol, e4, ctx);
  cov_mpoly_eval_fmpz(igusa_psi4(res), pol, M, ctx);
  cov_monomial(pol, e6, ctx);
  cov_mpoly_eval_fmpz(igusa_psi6(res), pol, M, ctx);
  cov_monomial(pol, e10, ctx);
  cov_mpoly_eval_fmpz(igusa_chi10(res), pol, M, ctx);
  cov_monomial(pol, e12, ctx);
  cov_mpoly_eval_fmpz(igusa_chi12(res), pol, M, ctx);

  if (z4) fmpz_zero(igusa_psi4(res));
  if (z6) fmpz_zero(igusa_psi6(res));
  if (z10) fmpz_zero(igusa_chi10(res));
  if (z12) fmpz_zero(igusa_chi12(res));

  _fmpz_vec_set(I, res, 4);
  cov_adjust_weights(weights, weights, I, 4);
  cov_normalize_fmpz(I, I, 4, weights);

  fmpz_mpoly_clear(pol, ctx);
  fmpz_mpoly_ctx_clear(ctx);
  flint_free(e4);
  flint_free(e6);
  flint_free(e10);
  flint_free(e12);
  _fmpz_vec_clear(res, 4);  
}


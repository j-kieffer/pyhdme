
#include "igusa.h"

void igusa_try_coordinate(fmpz_mpoly_t pol, slong wt, slong j,
			  const fmpz_mpoly_ctx_t ctx)
{
  slong k = igusa_nb_base_monomials(wt);
  slong bound;
  slong i;
  fmpz_t c;
  fmpz_mpoly_t term;

  fmpz_init(c);
  fmpz_mpoly_init(term, ctx);

  if (j < k)
    {
      igusa_base_monomial(pol, wt, j, ctx);
    }
  else
    {
      /* Decompose j in k numbers of roughly the same size */
      bound = FLINT_MAX(2, n_clog(j, k));
      fmpz_mpoly_zero(pol, ctx);
      for (i = 0; i < k; i++)
	{
	  fmpz_set_si(c, (j % bound) - (bound/2));
	  j = j/bound;
	  igusa_base_monomial(term, wt, i, ctx);
	  fmpz_mpoly_scalar_mul_fmpz(term, term, c, ctx);
	  fmpz_mpoly_add(pol, pol, term, ctx);
	}
    }

  fmpz_clear(c);
  fmpz_mpoly_clear(term, ctx);
}

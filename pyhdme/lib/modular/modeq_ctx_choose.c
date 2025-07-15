
#include "modular.h"

int modeq_ctx_choose(modeq_ctx_t ctx, acb_srcptr I, slong nb, slong prec)
{
  slong j, k, l;
  slong j0 = -1;
  acb_srcptr ptr; /* Do not initialize */
  slong weights[4] = IGUSA_HALFWEIGHTS;
  slong exps[4];
  fmpz_mpoly_t num, den;
  acb_ptr evden;
  acb_ptr evnum;
  int res;
  slong wt = 60;
  int v = get_modeq_verbose();

  fmpz_mpoly_init(num, modeq_ctx_ctx(ctx));
  fmpz_mpoly_init(den, modeq_ctx_ctx(ctx));
  evden = _acb_vec_init(nb);
  evnum = _acb_vec_init(nb);

  prec = modeq_ctx_prec(prec);

  /* Do we have: psi4 is never zero? */
  res = 1;
  for (j = 0; j < nb; j++)
    {
      if (acb_contains_zero(igusa_psi4(&I[4*j])))
	{
	  res = 0;
	  break;
	}
    }
  if (res) wt = 20;

  /* Do we have: psi4=psi6=0 and psi6=chi10=0 never happen? */
  res = 1;
  for (j = 0; j < nb; j++)
    {
      ptr = &I[4*j];
      if (acb_contains_zero(igusa_psi4(ptr))
	  && acb_contains_zero(igusa_psi6(ptr)))
	{
	  res = 0;
	  break;
	}
      if (acb_contains_zero(igusa_psi6(ptr))
	  && acb_contains_zero(igusa_chi10(ptr)))
	{
	  res = 0;
	  break;
	}
    }
  if (res && wt != 20) wt = 30;
  if (v) flint_printf("(modeq_ctx_choose) Take coordinates of weight %wd\n", wt);

  /* Set monomials in ctx */
  modeq_ctx_weight(ctx) = wt;
  modeq_ctx_nb(ctx) = igusa_nb_base_monomials(wt);
  for (k = 0; k < modeq_ctx_nb(ctx); k++)
    {
      igusa_base_exps(exps, wt, k);
      igusa_base_monomial(modeq_ctx_monomial(ctx, k), wt, k,
			  modeq_ctx_ctx(ctx));
    }

  /* Collect possibly equal pairs */
  for (j = 0; j < nb; j++)
    {
      j0 = -1;
      for (k = j+1; k < nb; k++)
	{
	  if (!cov_distinct(&I[4*j], &I[4*k], 4, weights, prec))
	    {
	      if (v && (j0==-1))
		{
		  flint_printf("(modeq_ctx_choose) Possibly isomorphic: %wd, %wd\n",
			       j, k);
		}
	      modeq_ctx_add_pair(ctx, j, k);
	      j0 = k;
	    }
	}
    }

  /* Choose denominator: should never vanish */
  for (j = 0; j < MODEQ_CTX_MAX_NB_COORDS; j++)
    {
      igusa_try_coordinate(den, wt, j, modeq_ctx_ctx(ctx));
      res = 1;
      for (k = 0; k < nb; k++)
	{
	  cov_mpoly_eval(&evden[k], den, &I[4*k], modeq_ctx_ctx(ctx), prec);
	  if (acb_contains_zero(&evden[k]))
	    {
	      res = 0;
	      break;
	    }
	}
      if (res)
	{
	  j0 = j;
	  break;
	}
    }
  
  if (res) fmpz_mpoly_set(modeq_ctx_den(ctx), den, modeq_ctx_ctx(ctx));
  if (res && v)
    {
      flint_printf("(modeq_ctx_choose) Denominator found: ");
      igusa_print_coordinate(modeq_ctx_den(ctx), modeq_ctx_ctx(ctx));
      flint_printf("\n");
    }
  else if (!res && v)
    {
      flint_printf("(modeq_ctx_choose) No suitable denominator found\n");
    }
  
  /* Choose numerator: num/den should be distinct on non-"equal" pairs */
  if (res)
    {
      for (j = 0; j < MODEQ_CTX_MAX_NB_COORDS; j++)
	{
	  if (j == j0) continue;
	  
	  igusa_try_coordinate(num, wt, j, modeq_ctx_ctx(ctx));	  
	  for (k = 0; k < nb; k++)
	    {
	      cov_mpoly_eval(&evnum[k], num, &I[4*k], modeq_ctx_ctx(ctx), prec);
	      acb_div(&evnum[k], &evnum[k], &evden[k], prec);
	    }
	  res = 1;
	  for (k = 0; k < nb; k++)
	    {
	      for (l = k+1; l < nb; l++)
		{
		  /* For performance, use that pairs are ordered? */
		  if (acb_overlaps(&evnum[k], &evnum[l])
		      && !modeq_ctx_is_pair(k, l, ctx))
		    {
		      res = 0;
		      break;
		    }
		}
	      if (!res) break;
	    }
	  if (res) break;
	}
      
      if (res) fmpz_mpoly_set(modeq_ctx_num(ctx), num, modeq_ctx_ctx(ctx));
      if (res && v)
	{
	  flint_printf("(modeq_ctx_choose) Numerator found: ");
	  igusa_print_coordinate(modeq_ctx_num(ctx), modeq_ctx_ctx(ctx));
	  flint_printf("\n");
	}
      else if (!res && v)
	{
	  flint_printf("(modeq_ctx_choose) No suitable numerator found\n");
	}
    }
    
  fmpz_mpoly_clear(num, modeq_ctx_ctx(ctx));
  fmpz_mpoly_clear(den, modeq_ctx_ctx(ctx));
  _acb_vec_clear(evden, nb);
  _acb_vec_clear(evnum, nb);
  return res;
}

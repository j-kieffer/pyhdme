
#include "modular.h"

int alt_2step_all_isog_Q(slong* nb_roots, fmpz* all_I, fmpz* I, slong ell)
{
  hecke_t H;
  modeq_t E;
  modeq_ctx_t ctx;
  slong nb_factors;
  slong max_nb_factors = siegel_nb_cosets(ell)/2;
  fmpz_poly_struct* factors;
  slong* mults;
  slong k;
  slong* indices;
  fmpz_mat_t L;
  slong add_nb;
  fmpz* add_I;
  int res, res2;
  int v = get_modeq_verbose();

  hecke_init(H, siegel_nb_cosets(ell));
  modeq_init(E);
  modeq_ctx_init(ctx);
  factors = flint_malloc(max_nb_factors * sizeof(fmpz_poly_struct));
  for (k = 0; k < max_nb_factors; k++) fmpz_poly_init(&factors[k]);
  mults = flint_malloc(max_nb_factors * sizeof(slong));
  indices = flint_malloc((ell+1) * sizeof(slong));
  fmpz_mat_init(L, 4, 4);
  add_I = _fmpz_vec_init(4 * siegel_nb_T1_cosets_with_line(ell));
  
  res = siegel_modeq_eval_with_hecke(E, ctx, H, I, ell);
  
  if (res)
    {
      *nb_roots = 0;
      alt_2step_factors(&nb_factors, factors, mults, E, ell);
      for (k = 0; k < nb_factors; k++)
	{
	  if (v) flint_printf("(alt_2step_all_isog_Q) Studying factor number %wd\n", k+1);
	  
	  /* Find out what the roots are; if failure, abort */
	  res = alt_2step_select_isog(indices, &factors[k], mults[k], H, ctx);
	  if (!res)
	    {
	      if (v) flint_printf("(alt_2step_all_isog_Q) Unable to isolate roots, abort.\n");
	      break;
	    }
	  
	  /* Find out what the stable line is; if failure, continue loop */
	  res2 = alt_2step_line(L, indices, fmpz_poly_degree(&factors[k]), H);
	  if (!res2)
	    {
	      if (v) flint_printf("(alt_2step_all_isog_Q) Found no stable line\n");	      
	      continue;
	    }
	  else if (v) flint_printf("(alt_2step_all_isog_Q) Stable line found; evaluating new modular equation\n");
	  
	  /* Evaluate new modular equation; if failure, abort */
	  res = alt_2step_modeq_with_line(E, ctx, L, I, ell);
	  if (!res) break;
	  
	  /* Compute and set roots */
	  modeq_all_isog_Q(&add_nb, add_I, E, ctx);
	  if (v) flint_printf("(alt_2step_all_isog_Q) Found %wd rational roots\n", add_nb);
	  
	  _fmpz_vec_set(&all_I[4*(*nb_roots)], add_I, 4*add_nb);
	  *nb_roots += add_nb;
	}
    }

  hecke_clear(H);
  modeq_clear(E);
  modeq_ctx_clear(ctx);
  for (k = 0; k < max_nb_factors; k++) fmpz_poly_clear(&factors[k]);
  flint_free(factors);
  flint_free(mults);
  flint_free(indices);
  fmpz_mat_clear(L);
  _fmpz_vec_clear(add_I, 4 * siegel_nb_T1_cosets_with_line(ell));
  
  return res;
}


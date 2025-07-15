
#include "hilbert.h"

int hilbert_linear_combination(fmpz* abcde, const acb_mat_t tau, slong delta, slong prec)
{
  fmpz_lll_t fl;
  fmpz_mat_t B, Bt;
  fmpz_mat_t U;
  fmpz_t discr;
  fmpz_t temp;
  acb_t coeff;
  slong exp = prec/2;
  slong k, l, m;
  int res;
  int verbose = get_hilbert_lll_verbose();

  fmpz_lll_context_init_default(fl);
  fmpz_mat_init(B, 7, 5);
  fmpz_mat_init(Bt, 5, 7);
  fmpz_mat_init(U, 5, 5);
  acb_init(coeff);
  fmpz_init(discr);
  fmpz_init(temp);

  /* Set up LLL matrix */
  fmpz_one(fmpz_mat_entry(B, 0, 0));
  fmpz_one(fmpz_mat_entry(B, 1, 1));
  fmpz_one(fmpz_mat_entry(B, 2, 2));
  fmpz_one(fmpz_mat_entry(B, 3, 3));
  fmpz_one(fmpz_mat_entry(B, 4, 4));

  acb_mul_2exp_si(coeff, acb_mat_entry(tau, 0, 0), exp);
  
  arf_get_fmpz(fmpz_mat_entry(B, 5, 0),
	       arb_midref(acb_realref(coeff)),
	       ARF_RND_NEAR);
  arf_get_fmpz(fmpz_mat_entry(B, 6, 0),
	       arb_midref(acb_imagref(coeff)),
	       ARF_RND_NEAR);
  
  acb_mul_2exp_si(coeff, acb_mat_entry(tau, 0, 1), exp);
  arf_get_fmpz(fmpz_mat_entry(B, 5, 1),
	       arb_midref(acb_realref(coeff)),
	       ARF_RND_NEAR);  
  arf_get_fmpz(fmpz_mat_entry(B, 6, 1),
	       arb_midref(acb_imagref(coeff)),
	       ARF_RND_NEAR);
  
  acb_mul_2exp_si(coeff, acb_mat_entry(tau, 1, 1), exp);
  arf_get_fmpz(fmpz_mat_entry(B, 5, 2),
	       arb_midref(acb_realref(coeff)),
	       ARF_RND_NEAR);  
  arf_get_fmpz(fmpz_mat_entry(B, 6, 2),
	       arb_midref(acb_imagref(coeff)),
	       ARF_RND_NEAR);

  acb_mat_det(coeff, tau, prec);
  acb_neg(coeff, coeff);
  acb_mul_2exp_si(coeff, coeff, exp);
  arf_get_fmpz(fmpz_mat_entry(B, 5, 3),
	       arb_midref(acb_realref(coeff)),
	       ARF_RND_NEAR);  
  arf_get_fmpz(fmpz_mat_entry(B, 6, 3),
	       arb_midref(acb_imagref(coeff)),
	       ARF_RND_NEAR);

  fmpz_one(fmpz_mat_entry(B, 5, 4));
  fmpz_mul_2exp(fmpz_mat_entry(B, 5, 4),
		fmpz_mat_entry(B, 5, 4),
		exp);

  /* Call LLL */
  if (verbose)
    {
      flint_printf("(hilbert_linear_combination) Matrix to reduce:\n");
      fmpz_mat_print_pretty(B);
      flint_printf("\n");
    }
  
  fmpz_mat_transpose(Bt, B);
  if (verbose)
    {
      flint_printf("(hilbert_linear_combination) Start LLL at precision %wd...", exp);
    }
  fflush(stdout);
  fmpz_lll(Bt, U, fl);
  if (verbose)
    {
      flint_printf(" done.\n"); fflush(stdout);
    }
  fmpz_mat_transpose(B, Bt);

  if (verbose)
    {
      flint_printf("(hilbert_linear_combination) First column should be small:\n");
      fmpz_mat_print_pretty(B);
      flint_printf("\n");
    }
  
  /* Get abcde from first column */
  for (k = 0; k < 5; k++)
    {
      fmpz_set(&abcde[k], fmpz_mat_entry(B, k, 0));
    }

  /* Check discriminant */
  fmpz_mul(discr, &abcde[1], &abcde[1]);
  fmpz_mul(temp, &abcde[0], &abcde[2]);
  fmpz_submul_ui(discr, temp, 4);
  fmpz_mul(temp, &abcde[3], &abcde[4]);
  fmpz_submul_ui(discr, temp, 4);
  res = fmpz_equal_si(discr, delta);

  if (!res)
    {
      /* Maybe a parasite relation, look at small linear combination
	 of two first columns having the right discriminant */
      for (l = -10; l <= 10; l++)
	{
	  for (m = -5; m <= 5; m++)
	    {
	      for (k = 0; k < 5; k++)
		{
		  fmpz_mul_si(&abcde[k], fmpz_mat_entry(B, k, 0), l);
      if (m > 0)
        fmpz_addmul_ui(&abcde[k], fmpz_mat_entry(B, k, 1), m);
      else
        fmpz_submul_ui(&abcde[k], fmpz_mat_entry(B, k, 1), -m);
		}
	      /* Check discriminant */
	      fmpz_mul(discr, &abcde[1], &abcde[1]);
	      fmpz_mul(temp, &abcde[0], &abcde[2]);
	      fmpz_submul_ui(discr, temp, 4);
	      fmpz_mul(temp, &abcde[3], &abcde[4]);
	      fmpz_submul_ui(discr, temp, 4);
	      res = fmpz_equal_si(discr, delta);
	      if (res) break;
	    }
	  if (res) break;
	}
    }
  
  fmpz_mat_clear(B);
  fmpz_mat_clear(Bt);
  fmpz_mat_clear(U);
  acb_clear(coeff);
  fmpz_clear(discr);
  fmpz_clear(temp);
  return res;
}

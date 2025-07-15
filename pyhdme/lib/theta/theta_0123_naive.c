
#include "theta.h"

/* Naive algorithm to compute theta_i for 0\leq i\leq 3.
   Notation is the same as Dupont's thesis, Alg. 15 (but B is not the same).
   
 */

int
theta_0123_naive(acb_ptr th, const acb_mat_t tau, slong prec)
{
  acb_t q1, q2, q3, q4;
  arb_t pi;
  
  fmpz_t B;
  slong B_slong;
  int res;
  
  acb_ptr Q1vec;
  acb_t a, b, a3, a4, b2, c2, s;
  acb_t Q2, Q3, Q4, S0, S1, S2, S3, A0, A1, A2, A3;
  slong m, n, pm1_n, pm1_m;
  int i;

  acb_init(q1);
  acb_init(q2);
  acb_init(q3);
  acb_init(q4);
  arb_init(pi);
  fmpz_init(B);

  acb_init(a);
  acb_init(b);
  acb_init(a3);
  acb_init(a4);
  acb_init(b2);
  acb_init(c2);
  acb_init(s);
  acb_init(Q2);
  acb_init(Q3);
  acb_init(Q4);
  acb_init(S0);
  acb_init(S1);
  acb_init(S2);
  acb_init(S3);
  acb_init(A0);
  acb_init(A1);
  acb_init(A2);
  acb_init(A3);

  /* Set q1, q2, q3, q4 */
  arb_const_pi(pi, prec);
  acb_mul_arb(q1, acb_mat_entry(tau, 0, 0), pi, prec);
  acb_mul_arb(q2, acb_mat_entry(tau, 1, 1), pi, prec);
  acb_mul_arb(q3, acb_mat_entry(tau, 0, 1), pi, prec);
  acb_mul_onei(q1, q1);
  acb_mul_onei(q2, q2);
  acb_mul_onei(q3, q3);
  acb_mul_2exp_si(q3, q3, 1);
  acb_neg(q4, q3);
  acb_exp(q1, q1, prec);
  acb_exp(q2, q2, prec);
  acb_exp(q3, q3, prec);
  acb_exp(q4, q4, prec);

  /* Adjust prec? */

  /* Compute B */
  res = theta_0123_naive_B(B, tau, prec);
  if (res) res = (fmpz_cmp_si(B, WORD_MAX) < 0); /* B+1 fits in a slong */
  if (res) B_slong = fmpz_get_si(B);

  if (res)
    {
      Q1vec = _acb_vec_init(B_slong + 2); /* 2 in case B=0 */
      
      /* Fill Q1vec */
      acb_set(a, q1);
      acb_sqr(b, q1, prec);
      acb_set(&Q1vec[1], q1);
      for (m = 2; m < B_slong+1; m++)
	{
	  acb_mul(a, a, b, prec);
	  acb_mul(&Q1vec[m], a, &Q1vec[m-1], prec);
	}
      
      /* Compute sums */
      /* Si are already zero */
      acb_one(Q2);
      acb_one(a3);
      acb_one(a4);
      acb_set(b2, q2);
      acb_sqr(c2, q2, prec);
      
      /* Outer loop */
      for (n = 1; n < B_slong+1; n++)
	{
	  acb_mul(a3, a3, q3, prec);
	  acb_mul(a4, a4, q4, prec);
	  acb_mul(Q2, Q2, b2, prec);
	  acb_mul(b2, b2, c2, prec);

	  pm1_n = (n%2 == 0 ? 1 : -1);
	  acb_add(S0, S0, &Q1vec[n], prec);
	  acb_addmul_si(S1, &Q1vec[n], pm1_n, prec);
	  acb_add(S2, S2, &Q1vec[n], prec);
	  acb_addmul_si(S3, &Q1vec[n], pm1_n, prec);

	  acb_one(Q3);
	  acb_one(Q4);
	  acb_one(A0);
	  acb_one(A1);
	  acb_set_si(A2, pm1_n);
	  acb_set_si(A3, pm1_n);
	  
	  /* Inner loop */
	  for (m = 1; m < B_slong+1; m++)
	    {
	      acb_mul(Q3, Q3, a3, prec);
	      acb_mul(Q4, Q4, a4, prec);
	      acb_add(s, Q3, Q4, prec);
	      acb_mul(s, s, &Q1vec[m], prec);

	      pm1_m = (m%2 == 0 ? 1 : -1);
	      acb_add(A0, A0, s, prec);
	      acb_addmul_si(A1, s, pm1_m, prec);
	      acb_addmul_si(A2, s, pm1_n, prec);
	      acb_addmul_si(A3, s, pm1_n * pm1_m, prec);
	    } /* End inner loop */
	  acb_addmul(S0, Q2, A0, prec);
	  acb_addmul(S1, Q2, A1, prec);
	  acb_addmul(S2, Q2, A2, prec);
	  acb_addmul(S3, Q2, A3, prec);
	} /* End outer loop */
      
      /* Clear Q1vec */
      _acb_vec_clear(Q1vec, B_slong + 2);

      /* Return 1 + 2*Si */
      /* Following our conventions, 1 and 2 are inverted. */
      acb_one(&th[0]);
      acb_addmul_si(&th[0], S0, 2, prec);
      acb_one(&th[2]);
      acb_addmul_si(&th[2], S1, 2, prec);
      acb_one(&th[1]);
      acb_addmul_si(&th[1], S2, 2, prec);
      acb_one(&th[3]);
      acb_addmul_si(&th[3], S3, 2, prec);

      /* Add error from tail of series */
      for (i = 0; i < 4; i++)
	{
	  arb_add_error_2exp_si(acb_realref(&th[i]), -prec);
	  arb_add_error_2exp_si(acb_imagref(&th[i]), -prec);
	}
    }
  
  acb_clear(q1);
  acb_clear(q2);
  acb_clear(q3);
  acb_clear(q4);
  arb_clear(pi);
  fmpz_clear(B);
  
  acb_clear(a);
  acb_clear(b);
  acb_clear(a3);
  acb_clear(a4);
  acb_clear(b2);
  acb_clear(c2);
  acb_clear(s);
  acb_clear(Q2);
  acb_clear(Q3);
  acb_clear(Q4);
  acb_clear(S0);
  acb_clear(S1);
  acb_clear(S2);
  acb_clear(S3);
  acb_clear(A0);
  acb_clear(A1);
  acb_clear(A2);
  acb_clear(A3);
  
  return res;
}

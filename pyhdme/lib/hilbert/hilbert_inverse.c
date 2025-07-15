
#include "hilbert.h"

static void
hilbert_bw_M0(fmpz_mat_t m, fmpz* abcde)
{
  fmpz_t g0, alpha, beta, one;
  fmpz_mat_t M0, R0, R1;
  
  fmpz_init(g0);
  fmpz_init(alpha);
  fmpz_init(beta);
  fmpz_init(one);
  fmpz_mat_init(M0, 4, 4);
  fmpz_mat_init(R0, 4, 4);
  fmpz_mat_init(R1, 4, 4);

  fmpz_xgcd(g0, alpha, beta, &abcde[4], &abcde[2]);
  fmpz_neg(beta, beta);

  fmpz_mat_zero(M0);
  fmpz_one(fmpz_mat_entry(M0, 0, 0));
  fmpz_one(fmpz_mat_entry(M0, 0, 2));
  fmpz_set(fmpz_mat_entry(M0, 1, 1), alpha);
  fmpz_set(fmpz_mat_entry(M0, 1, 3), beta);
  fmpz_set_si(fmpz_mat_entry(M0, 2, 0), -1);
  fmpz_divexact(fmpz_mat_entry(M0, 3, 1), &abcde[2], g0);
  fmpz_divexact(fmpz_mat_entry(M0, 3, 3), &abcde[4], g0);

  fmpz_mat_set(m, M0);
  if (!fmpz_mat_is_symplectic(m))
    {
      flint_printf("(hilbert_bw_M0) Incorrect m\n");
      fmpz_mat_print(m);
      fflush(stdout);
      flint_abort();
    }

  fmpz_mat_zero(R0);
  fmpz_set(fmpz_mat_entry(R0, 0, 1), &abcde[0]);
  fmpz_set(fmpz_mat_entry(R0, 0, 3), &abcde[3]);
  fmpz_neg(fmpz_mat_entry(R0, 1, 0), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R0, 1, 1), &abcde[1]);
  fmpz_neg(fmpz_mat_entry(R0, 1, 2), &abcde[3]);
  fmpz_set(fmpz_mat_entry(R0, 2, 1), &abcde[4]);
  fmpz_neg(fmpz_mat_entry(R0, 2, 3), &abcde[2]);
  fmpz_neg(fmpz_mat_entry(R0, 3, 0), &abcde[4]);
  fmpz_set(fmpz_mat_entry(R0, 3, 2), &abcde[0]);
  fmpz_set(fmpz_mat_entry(R0, 3, 3), &abcde[1]);

  fmpz_mat_transpose(M0, M0);
  fmpz_mat_mul(R0, R0, M0);
  fmpz_mat_inv(M0, one, M0);
  fmpz_mat_mul(R0, M0, R0);
  
  if (!fmpz_is_one(one))
    {
      flint_printf("(hilbert_bw_M0) Incorrect inversion\n");
      fflush(stdout);
      flint_abort();
    }
  
  fmpz_set(&abcde[0], fmpz_mat_entry(R0, 0, 1));
  fmpz_set(&abcde[1], fmpz_mat_entry(R0, 1, 1));
  fmpz_neg(&abcde[2], fmpz_mat_entry(R0, 1, 0));
  fmpz_zero(&abcde[3]);
  fmpz_set(&abcde[4], fmpz_mat_entry(R0, 2, 1));

  /* Check that new R0 is indeed the right R1 */
  fmpz_mat_zero(R1);
  fmpz_set(fmpz_mat_entry(R1, 0, 1), &abcde[0]);
  fmpz_neg(fmpz_mat_entry(R1, 1, 0), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R1, 1, 1), &abcde[1]);
  fmpz_set(fmpz_mat_entry(R1, 2, 1), &abcde[4]);
  fmpz_neg(fmpz_mat_entry(R1, 2, 3), &abcde[2]);
  fmpz_neg(fmpz_mat_entry(R1, 3, 0), &abcde[4]);
  fmpz_set(fmpz_mat_entry(R1, 3, 2), &abcde[0]);
  fmpz_set(fmpz_mat_entry(R1, 3, 3), &abcde[1]);
  if (!fmpz_mat_equal(R0, R1))
    {
      flint_printf("(hilbert_bw_M0) Incorrect R1\n");
      fmpz_mat_print_pretty(R0); flint_printf("\n");
      fmpz_mat_print_pretty(R1); flint_printf("\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_clear(g0);
  fmpz_clear(alpha);
  fmpz_clear(beta);
  fmpz_clear(one);
  fmpz_mat_clear(M0);
  fmpz_mat_clear(R0);
  fmpz_mat_clear(R1);
}


/* Change from Birkenhake-Wilhelm: we only need p coprime to c1. */
static void
hilbert_bw_M1(fmpz_mat_t m, fmpz* abcde)
{
  fmpz_t g1, n, p, absc1, pstep;
  fmpz_mat_t M1;
  
  fmpz_init(g1);
  fmpz_init(n);
  fmpz_init(p);
  fmpz_init(absc1);
  fmpz_init(pstep);
  fmpz_mat_init(M1, 4, 4);

  fmpz_gcd(g1, &abcde[0], &abcde[4]);

  fmpz_divexact(p, &abcde[4], g1);
  fmpz_divexact(pstep, &abcde[0], g1);
  fmpz_abs(pstep, pstep);
  fmpz_abs(absc1, &abcde[2]);
  
  /* Use g1 in coprimality test */
  fmpz_zero(n);
  fmpz_gcd(g1, p, absc1);
  while (!fmpz_is_one(g1))
    {
      fmpz_add_si(n, n, 1);
      fmpz_add(p, p, pstep);
      fmpz_gcd(g1, p, absc1);
    }
  
  fmpz_mat_one(M1);
  fmpz_neg(fmpz_mat_entry(M1, 0, 2), n);
  fmpz_mat_set(m, M1);
  if (!fmpz_mat_is_symplectic(m))
    {
      flint_printf("(hilbert_bw_M1) Incorrect m\n");
      fmpz_mat_print(m);
      fflush(stdout);
      flint_abort();
    }

  fmpz_mul(n, n, &abcde[0]);
  fmpz_add(&abcde[4], &abcde[4], n);

  fmpz_clear(g1);
  fmpz_clear(n);
  fmpz_clear(p);
  fmpz_clear(absc1);
  fmpz_clear(pstep);
  fmpz_mat_clear(M1);
}

static void
hilbert_bw_M2(fmpz_mat_t m, fmpz* abcde)
{
  fmpz_t g2, gamma, delta, one;
  fmpz_mat_t M2, R2, R3;
  
  fmpz_init(g2);
  fmpz_init(gamma);
  fmpz_init(delta);
  fmpz_init(one);
  fmpz_mat_init(M2, 4, 4);
  fmpz_mat_init(R2, 4, 4);
  fmpz_mat_init(R3, 4, 4);

  fmpz_xgcd(g2, gamma, delta, &abcde[4], &abcde[2]);
  fmpz_neg(delta, delta);

  fmpz_mat_zero(M2);
  fmpz_mul(fmpz_mat_entry(M2, 0, 0), &abcde[0], gamma);
  fmpz_divexact(fmpz_mat_entry(M2, 0, 0), fmpz_mat_entry(M2, 0, 0), g2);
  fmpz_one(fmpz_mat_entry(M2, 0, 2));
  fmpz_set(fmpz_mat_entry(M2, 1, 1), gamma);
  fmpz_set(fmpz_mat_entry(M2, 1, 3), delta);
  fmpz_set_si(fmpz_mat_entry(M2, 2, 0), -1);
  fmpz_divexact(fmpz_mat_entry(M2, 3, 1), &abcde[2], g2);
  fmpz_divexact(fmpz_mat_entry(M2, 3, 3), &abcde[4], g2);

  fmpz_mat_set(m, M2);
  if (!fmpz_mat_is_symplectic(m))
    {
      flint_printf("(hilbert_bw_M2) Incorrect m\n");
      fmpz_mat_print(m);
      fflush(stdout);
      flint_abort();
    }

  fmpz_mat_zero(R2);
  fmpz_set(fmpz_mat_entry(R2, 0, 1), &abcde[0]);
  fmpz_neg(fmpz_mat_entry(R2, 1, 0), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R2, 1, 1), &abcde[1]);
  fmpz_set(fmpz_mat_entry(R2, 2, 1), &abcde[4]);
  fmpz_neg(fmpz_mat_entry(R2, 2, 3), &abcde[2]);
  fmpz_neg(fmpz_mat_entry(R2, 3, 0), &abcde[4]);
  fmpz_set(fmpz_mat_entry(R2, 3, 2), &abcde[0]);
  fmpz_set(fmpz_mat_entry(R2, 3, 3), &abcde[1]);

  fmpz_mat_transpose(M2, M2);
  fmpz_mat_mul(R2, R2, M2);
  fmpz_mat_inv(M2, one, M2);
  fmpz_mat_mul(R2, M2, R2);
  
  if (!fmpz_is_one(one))
    {
      flint_printf("(hilbert_bw_M2) Incorrect inversion\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_set(&abcde[0], fmpz_mat_entry(R2, 0, 1));
  fmpz_set(&abcde[1], fmpz_mat_entry(R2, 1, 1));
  fmpz_neg(&abcde[2], fmpz_mat_entry(R2, 1, 0));
  fmpz_zero(&abcde[4]);
  
  /* Check that new R2 is indeed the right R3 */
  fmpz_mat_zero(R3);
  fmpz_set(fmpz_mat_entry(R3, 0, 1), &abcde[0]);
  fmpz_neg(fmpz_mat_entry(R3, 1, 0), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R3, 1, 1), &abcde[1]);
  fmpz_neg(fmpz_mat_entry(R3, 2, 3), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R3, 3, 2), &abcde[0]);
  fmpz_set(fmpz_mat_entry(R3, 3, 3), &abcde[1]);
  if (!fmpz_mat_equal(R2, R3))
    {
      flint_printf("(hilbert_bw_M2) Incorrect R3\n");
      fmpz_mat_print_pretty(R2); flint_printf("\n");
      fmpz_mat_print_pretty(R3); flint_printf("\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_clear(g2);
  fmpz_clear(gamma);
  fmpz_clear(delta);
  fmpz_clear(one);
  fmpz_mat_clear(M2);
  fmpz_mat_clear(R2);
  fmpz_mat_clear(R3);
}

static void
hilbert_bw_M3(fmpz_mat_t m, fmpz* abcde)
{
  fmpz_t g3, epsilon, eta, one;
  fmpz_mat_t M3, R3, R4;

  fmpz_init(g3);
  fmpz_init(epsilon);
  fmpz_init(eta);
  fmpz_init(one);
  fmpz_mat_init(M3, 4, 4);
  fmpz_mat_init(R3, 4, 4);
  fmpz_mat_init(R4, 4, 4);

  fmpz_xgcd(g3, epsilon, eta, &abcde[0], &abcde[2]);

  fmpz_mat_zero(M3);
  fmpz_one(fmpz_mat_entry(M3, 0, 0));
  fmpz_one(fmpz_mat_entry(M3, 0, 2));
  fmpz_divexact(fmpz_mat_entry(M3, 1, 1), &abcde[2], g3);
  fmpz_neg(fmpz_mat_entry(M3, 1, 3), &abcde[0]);
  fmpz_divexact(fmpz_mat_entry(M3, 1, 3), fmpz_mat_entry(M3, 1, 3), g3);
  fmpz_mul(fmpz_mat_entry(M3, 2, 0), epsilon, &abcde[0]);
  fmpz_neg(fmpz_mat_entry(M3, 2, 0), fmpz_mat_entry(M3, 2, 0));
  fmpz_divexact(fmpz_mat_entry(M3, 2, 0), fmpz_mat_entry(M3, 2, 0), g3);
  fmpz_mul(fmpz_mat_entry(M3, 2, 2), eta, &abcde[2]);
  fmpz_divexact(fmpz_mat_entry(M3, 2, 2), fmpz_mat_entry(M3, 2, 2), g3);
  fmpz_set(fmpz_mat_entry(M3, 3, 1), epsilon);
  fmpz_set(fmpz_mat_entry(M3, 3, 3), eta);

  fmpz_mat_set(m, M3);
  if (!fmpz_mat_is_symplectic(m))
    {
      flint_printf("(hilbert_bw_M3) Incorrect m\n");
      fmpz_mat_print(m);
      fflush(stdout);
      flint_abort();
    }

  fmpz_mat_zero(R3);
  fmpz_set(fmpz_mat_entry(R3, 0, 1), &abcde[0]);
  fmpz_neg(fmpz_mat_entry(R3, 1, 0), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R3, 1, 1), &abcde[1]);
  fmpz_neg(fmpz_mat_entry(R3, 2, 3), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R3, 3, 2), &abcde[0]);
  fmpz_set(fmpz_mat_entry(R3, 3, 3), &abcde[1]);

  fmpz_mat_transpose(M3, M3);
  fmpz_mat_mul(R3, R3, M3);
  fmpz_mat_inv(M3, one, M3);
  fmpz_mat_mul(R3, M3, R3);
  
  if (!fmpz_is_one(one))
    {
      flint_printf("(hilbert_bw_M3) Incorrect inversion\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_neg(&abcde[2], fmpz_mat_entry(R3, 1, 0));
  fmpz_divexact(&abcde[0], fmpz_mat_entry(R3, 0, 1), &abcde[2]);
  fmpz_set(&abcde[1], fmpz_mat_entry(R3, 1, 1));

  /* Check that new R3 is indeed the right R4 */
  fmpz_mat_zero(R4);
  fmpz_mul(fmpz_mat_entry(R4, 0, 1), &abcde[0], &abcde[2]);
  fmpz_neg(fmpz_mat_entry(R4, 1, 0), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R4, 1, 1), &abcde[1]);
  fmpz_neg(fmpz_mat_entry(R4, 2, 3), &abcde[2]);
  fmpz_mul(fmpz_mat_entry(R4, 3, 2), &abcde[0], &abcde[2]);
  fmpz_set(fmpz_mat_entry(R4, 3, 3), &abcde[1]);
  if (!fmpz_mat_equal(R3, R4))
    {
      flint_printf("(hilbert_bw_M3) Incorrect R4\n");
      fmpz_mat_print_pretty(R3); flint_printf("\n");
      fmpz_mat_print_pretty(R4); flint_printf("\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_clear(g3);
  fmpz_clear(epsilon);
  fmpz_clear(eta);
  fmpz_clear(one);
  fmpz_mat_clear(M3);
  fmpz_mat_clear(R3);
  fmpz_mat_clear(R4);
}

static void
hilbert_bw_M4(fmpz_mat_t m, fmpz* abcde)
{
  fmpz_t lhs, mu, nu, one;
  fmpz_mat_t M4, M4prime, R4, R5;

  fmpz_init(lhs);
  fmpz_init(mu);
  fmpz_init(nu);
  fmpz_init(one);
  fmpz_mat_init(M4, 4, 4);
  fmpz_mat_init(M4prime, 4, 4);
  fmpz_mat_init(R4, 4, 4);
  fmpz_mat_init(R5, 4, 4);

  fmpz_add_si(lhs, &abcde[0], 1);
  fmpz_mul(lhs, lhs, &abcde[2]);
  fmpz_add(lhs, lhs, &abcde[1]);
  fmpz_xgcd(one, mu, nu, lhs, &abcde[2]);

  if (!fmpz_is_one(one))
    {
      flint_printf("(hilbert_bw_M3) gcd is not one:\n");
      fmpz_print(one); flint_printf("\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_mat_zero(M4);
  fmpz_one(fmpz_mat_entry(M4, 0, 0));
  fmpz_one(fmpz_mat_entry(M4, 0, 2));
  fmpz_set(fmpz_mat_entry(M4, 1, 1), &abcde[2]);
  fmpz_neg(fmpz_mat_entry(M4, 1, 3), lhs);
  fmpz_mul(fmpz_mat_entry(M4, 2, 0), mu, lhs);
  fmpz_neg(fmpz_mat_entry(M4, 2, 0), fmpz_mat_entry(M4, 2, 0));
  fmpz_mul(fmpz_mat_entry(M4, 2, 2), nu, &abcde[2]);
  fmpz_set(fmpz_mat_entry(M4, 3, 1), mu);
  fmpz_set(fmpz_mat_entry(M4, 3, 3), nu);

  fmpz_mat_set(m, M4);
  if (!fmpz_mat_is_symplectic(m))
    {
      flint_printf("(hilbert_bw_M4) Incorrect m\n");
      fmpz_mat_print(m);
      fflush(stdout);
      flint_abort();
    }

  fmpz_mat_one(M4prime);
  fmpz_set_si(fmpz_mat_entry(M4prime, 1, 0), -1);
  fmpz_one(fmpz_mat_entry(M4prime, 2, 3));
  fmpz_mat_mul(M4, M4, M4prime);
  
  fmpz_mat_set(m, M4);
  if (!fmpz_mat_is_symplectic(m))
    {
      flint_printf("(hilbert_bw_M4) Incorrect m\n");
      fmpz_mat_print(m);
      fflush(stdout);
      flint_abort();
    }

  fmpz_mat_zero(R4);
  fmpz_mul(fmpz_mat_entry(R4, 0, 1), &abcde[0], &abcde[2]);
  fmpz_neg(fmpz_mat_entry(R4, 1, 0), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R4, 1, 1), &abcde[1]);
  fmpz_neg(fmpz_mat_entry(R4, 2, 3), &abcde[2]);
  fmpz_mul(fmpz_mat_entry(R4, 3, 2), &abcde[0], &abcde[2]);
  fmpz_set(fmpz_mat_entry(R4, 3, 3), &abcde[1]);

  fmpz_mat_transpose(M4, M4);
  fmpz_mat_mul(R4, R4, M4);
  fmpz_mat_inv(M4, one, M4);
  fmpz_mat_mul(R4, M4, R4);
  fmpz_mat_one(M4prime);
  fmpz_mat_scalar_mul_fmpz(M4prime, M4prime, &abcde[2]);
  fmpz_mat_add(R4, R4, M4prime);

  if (!fmpz_is_one(one))
    {
      flint_printf("(hilbert_bw_M1) Incorrect inversion\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_set(&abcde[0], fmpz_mat_entry(R4, 0, 1));
  fmpz_set(&abcde[1], fmpz_mat_entry(R4, 1, 1));
  fmpz_one(&abcde[2]);
  
  /* Check that new R4 is indeed the right R5 */
  fmpz_mat_zero(R5);
  fmpz_set(fmpz_mat_entry(R5, 0, 1), &abcde[0]);
  fmpz_neg(fmpz_mat_entry(R5, 1, 0), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R5, 1, 1), &abcde[1]);
  fmpz_neg(fmpz_mat_entry(R5, 2, 3), &abcde[2]);
  fmpz_set(fmpz_mat_entry(R5, 3, 2), &abcde[0]);
  fmpz_set(fmpz_mat_entry(R5, 3, 3), &abcde[1]);
  if (!fmpz_mat_equal(R4, R5))
    {
      flint_printf("(hilbert_bw_M4) Incorrect R5\n");
      fmpz_mat_print_pretty(R4); flint_printf("\n");
      fmpz_mat_print_pretty(R5); flint_printf("\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_clear(lhs);
  fmpz_clear(mu);
  fmpz_clear(nu);
  fmpz_clear(one);
  fmpz_mat_clear(M4);
  fmpz_mat_clear(M4prime);
  fmpz_mat_clear(R4);
  fmpz_mat_clear(R5);
}

/* Change from Birkenhake-Wilhelm: we want b'=-1 if b5 is odd. */
static void
hilbert_bw_M5(fmpz_mat_t m, fmpz* abcde)
{
  fmpz_t t;
  fmpz_mat_t M5;

  fmpz_init(t);
  fmpz_mat_init(M5, 4, 4);
  
  if (fmpz_is_even(&abcde[1]))
    {
      fmpz_divexact_si(t, &abcde[1], 2);
    }
  else
    {
      fmpz_add_si(t, &abcde[1], 1);
      fmpz_divexact_si(t, t, 2);
    }

  fmpz_mat_one(M5);
  fmpz_set(fmpz_mat_entry(M5, 1, 0), t);
  fmpz_neg(fmpz_mat_entry(M5, 2, 3), t);

  fmpz_mat_set(m, M5);
  if (!fmpz_mat_is_symplectic(m))
    {
      flint_printf("(hilbert_bw_M5) Incorrect m\n");
      fmpz_mat_print(m);
      fflush(stdout);
      flint_abort();
    }

  fmpz_clear(t);
  fmpz_mat_clear(M5);
}

int hilbert_inverse(acb_ptr t, fmpz_mat_t eta, const acb_mat_t tau,
		    slong delta, slong prec)
{
  fmpz* abcde;
  int res;
  fmpz_mat_t m;
  acb_mat_t im;
  acb_mat_t R;

  abcde = _fmpz_vec_init(5);
  fmpz_mat_init(m, 4, 4);
  acb_mat_init(im, 2, 2);
  acb_mat_init(R, 2, 2);

  hilbert_R(R, delta, prec);
  res = hilbert_linear_combination(abcde, tau, delta, prec);
  if (res)
    {
      hilbert_bw_M0(m, abcde);
      fmpz_mat_set(eta, m);
      hilbert_bw_M1(m, abcde);
      fmpz_mat_mul(eta, m, eta);
      hilbert_bw_M2(m, abcde);
      fmpz_mat_mul(eta, m, eta);
      hilbert_bw_M3(m, abcde);
      fmpz_mat_mul(eta, m, eta);
      hilbert_bw_M4(m, abcde);
      fmpz_mat_mul(eta, m, eta);
      hilbert_bw_M5(m, abcde);
      fmpz_mat_mul(eta, m, eta);

      siegel_transform(im, eta, tau, prec);
      acb_mat_inv(R, R, prec);
      acb_mat_mul(im, im, R, prec);
      acb_mat_transpose(R, R);
      acb_mat_mul(im, R, im, prec);
      
      if (!acb_contains_zero(acb_mat_entry(im, 0, 1)))
	{
	  flint_printf("(hilbert_inverse) Warning: could not bring tau to normal form\n");
	  acb_mat_printd(tau, 30); flint_printf("\n");
	  res = 0;
	}
    }
  if (res)
    {
      acb_set(&t[0], acb_mat_entry(im, 0, 0));
      acb_set(&t[1], acb_mat_entry(im, 1, 1));
    }  
  
  _fmpz_vec_clear(abcde, 5);
  fmpz_mat_clear(m);
  acb_mat_clear(im);
  acb_mat_clear(R);
  return res;
}

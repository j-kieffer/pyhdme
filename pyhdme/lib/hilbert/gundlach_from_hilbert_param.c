
#include "hilbert.h"

void gundlach_from_hilbert_param(fmpz* G, fmpq* mn, slong delta)
{
  fmpq_t g;
  fmpq_t h;
  fmpq_t temp, temp2;
  slong weights[3] = {1, 3, 5};

  if (delta != 5)
    {
      flint_printf("(gundlach_from_hilbert_param) Gundlach invariants only implemented for delta = 5\n");
      fflush(stdout);
      flint_abort();
    }

  fmpq_init(g);
  fmpq_init(h);
  fmpq_init(temp);
  fmpq_init(temp2);

  fmpq_pow_si(g, &mn[0], 2);
  fmpq_pow_si(temp, &mn[1], 2);
  fmpq_mul_si(temp, temp, -5);
  fmpq_add(g, g, temp);
  fmpq_add_si(g, g, -9);
  fmpq_set_si(temp, 30, 1);
  fmpq_div(g, g, temp);

  fmpq_mul_si(h, &mn[0], 3);
  fmpq_mul_si(temp, g, 10);
  fmpq_add_si(temp, temp, 3);
  fmpq_mul(h, h, temp);
  fmpq_mul_si(temp, g, 15);
  fmpq_add_si(temp, temp, 2);
  fmpq_mul(h, h, temp);
  fmpq_pow_si(temp, g, 2);
  fmpq_mul_si(temp, temp, 250);
  fmpq_mul_si(temp2, g, 75);
  fmpq_add(temp, temp, temp2);
  fmpq_add_si(temp, temp, 6);
  fmpq_mul_si(temp, temp, 9);
  fmpq_add(h, h, temp);
  fmpq_set_si(temp, 6250, 1);
  fmpq_div(h, h, temp);

  /* Lift to Gundlach covariants */
  /* See Milio--Robert, "Modular polynomials on Hilbert surfaces" §A.2 and §A.4 */
  /* Can take G2 = -6g, F6 = h, F10 = h^2; rescale to ensure they are integers */
  fmpz_mul_si(&G[0], fmpq_numref(g), -6);
  fmpz_set(&G[1], fmpq_numref(h));
  fmpz_pow_ui(&G[2], &G[1], 2);

  cov_rescale_fmpz(G, G, fmpq_denref(g), 3, weights);
  fmpz_divexact(&G[0], &G[0], fmpq_denref(g));
  cov_rescale_fmpz(G, G, fmpq_denref(h), 3, weights);
  fmpz_divexact(&G[1], &G[1], fmpq_denref(h));
  fmpz_divexact(&G[2], &G[2], fmpq_denref(h));
  fmpz_divexact(&G[2], &G[2], fmpq_denref(h));
  cov_normalize_fmpz(G, G, 3, weights);

  fmpq_clear(g);
  fmpq_clear(h);
  fmpq_clear(temp);
  fmpq_clear(temp2);
}

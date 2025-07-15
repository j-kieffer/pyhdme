
#include "modular.h"

void modeq_scalar(acb_t c, const hecke_t H, fmpz* I,
		  const modeq_ctx_t ctx, slong prec)
{
  slong wt = modeq_ctx_weight(ctx);
  slong weights[4] = IGUSA_WEIGHTS;
  acb_t s;
  fmpq_t scal;
  slong k;

  acb_init(s);
  fmpq_init(scal);
  
  /* Step 1: product of all stardets to the correct weight, to get
     a Siegel modular form */
  acb_one(c);
  for (k = 0; k < hecke_nb(H); k++)
    {
      acb_inv(s, hecke_stardet(H, k), prec);
      acb_mul(c, c, s, prec);
    }
  
  /* Step 2: normalize Hecke operator */
  acb_set_fmpz(s, hecke_norm_all(H));
  acb_mul(c, c, s, prec);
  
  /* Step 3: multiply by scaling factor between hecke_I_tau and I */
  cov_find_rescaling(s, hecke_I_tau(H), I, 4, weights, prec);
  acb_pow_si(s, s, -hecke_nb(H), prec);
  acb_mul(c, c, s, prec);

  /* Step 4: compute power according to weight */
  acb_pow_si(c, c, wt, prec);

  /* Step 5: additional factor according to additional MF generators */
  igusa_make_integral(scal, I, hecke_nb(H) * wt, FLINT_MAX(hecke_prod_ec(H) - 1, 0) * wt);
  acb_set_fmpq(s, scal, prec);
  acb_mul(c, c, s, prec);

  /* Step 6: divide by power of chi10 if appropriate */
  if (modeq_ctx_has_chi10_factor(ctx) && hecke_prod_ec(H) >= 1)
    {
      acb_set_fmpz(s, igusa_chi10(I));
      acb_pow_si(s, s, 1 - hecke_prod_ec(H), prec);
      acb_mul(c, c, s, prec);
    }

  acb_clear(s);
  fmpq_clear(scal);
}

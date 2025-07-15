
#include "igusa.h"

/* Following Bolza, "On binary sextics with linear transformations
   into themselves", p. 51 ff. */   

static void fmpq_div_si(fmpq_t r, const fmpq_t n, slong d)
{
  fmpz_t c;
  fmpz_init(c);
  fmpz_set_si(c, d);
  fmpq_div_fmpz(r, n, c);
  fmpz_clear(c);
}

static int bolza_conditions_12(fmpq* ABCD)
{
  return (fmpq_is_zero(&ABCD[1])
	  && fmpq_is_zero(&ABCD[2])
	  && fmpq_is_zero(&ABCD[3]));
}

static int bolza_conditions_11(fmpq* ABCD)
{
  fmpq_t t1, t2, temp;
  int res;
  
  fmpq_init(t1);
  fmpq_init(t2);
  fmpq_init(temp);
  
  fmpq_mul(temp, &ABCD[0], &ABCD[0]);
  fmpq_mul_si(t1, &ABCD[1], 6);
  fmpq_sub(t1, t1, temp);

  fmpq_mul(temp, &ABCD[0], &ABCD[1]);
  fmpq_div_si(temp, temp, 3);
  fmpq_mul_si(t2, &ABCD[2], 2);
  fmpq_add(t2, t2, temp);

  res = (fmpq_is_zero(t1)
	 && fmpq_is_zero(t2)
	 && fmpq_is_zero(&ABCD[3]));

  fmpq_clear(t1);
  fmpq_clear(t2);
  fmpq_clear(temp);
  return res;
}

static int bolza_conditions_19(acb_t a, fmpq* ABCD, slong prec)
{
  fmpq_t t1, t2, t3, t4, temp;
  int res;

  fmpq_init(t1);
  fmpq_init(t2);
  fmpq_init(t3);
  fmpq_init(t4);
  fmpq_init(temp);

  fmpq_pow_si(t1, &ABCD[1], 3);
  fmpq_div_si(t1, t1, -6);
  fmpq_pow_si(temp, &ABCD[2], 2);
  fmpq_add(t1, t1, temp);
  
  fmpq_mul(t4, &ABCD[0], &ABCD[1]);
  fmpq_mul_si(temp, &ABCD[2], 6);
  fmpq_add(t4, t4, temp);

  fmpq_mul(t2, t4, &ABCD[1]);
  fmpq_mul_si(t2, t2, -2);
  fmpq_mul_si(temp, &ABCD[3], 9);
  fmpq_add(t2, t2, temp);

  fmpq_mul_si(t3, &ABCD[2], -15);
  fmpq_mul(temp, &ABCD[0], &ABCD[1]);
  fmpq_mul_si(temp, temp, 2);
  fmpq_add(t3, t3, temp);

  res = (fmpq_is_zero(t1) && fmpq_is_zero(t2));
  if (res)
    {
      if (fmpq_is_zero(&ABCD[3]) || fmpq_is_zero(t3))
	{
	  flint_printf("(mestre_fmpz: bolza_conditions_19) Error: unexpected vanishing\n");
	  fflush(stdout);
	  flint_abort();
	}
      fmpq_div(t4, t4, t3);
      fmpq_neg(t4, t4);
      acb_set_fmpq(a, t4, prec);
      borchardt_sqrt(a, a, prec);
      acb_mul_si(a, a, 10, prec);
    }

  fmpq_clear(t1);
  fmpq_clear(t2);
  fmpq_clear(t3);
  fmpq_clear(t4);
  fmpq_clear(temp);
  return res;
}

static int bolza_conditions_23(acb_t a, fmpq* ABCD, slong prec)
{
  fmpq_t t1, t2, t3, t4, temp;
  int res;

  fmpq_init(t1);
  fmpq_init(t2);
  fmpq_init(t3);
  fmpq_init(t4);
  fmpq_init(temp);

  fmpq_pow_si(t1, &ABCD[1], 2);
  fmpq_mul_si(temp, &ABCD[0], 3);
  fmpq_mul(t1, t1, temp);
  fmpq_mul(temp, &ABCD[1], &ABCD[2]);
  fmpq_mul_si(temp, temp, -6);
  fmpq_add(t1, t1, temp);
  fmpq_pow_si(temp, &ABCD[0], 2);
  fmpq_mul(temp, temp, &ABCD[2]);
  fmpq_mul_si(temp, temp, 4);
  fmpq_add(t1, t1, temp);
  fmpq_mul_si(temp, &ABCD[3], -18);
  fmpq_add(t1, t1, temp);

  fmpq_pow_si(t2, &ABCD[1], 3);
  fmpq_mul_si(t2, t2, 4);
  fmpq_mul(temp, &ABCD[0], &ABCD[1]);
  fmpq_mul(temp, temp, &ABCD[2]);
  fmpq_mul_si(temp, temp, 5);
  fmpq_add(t2, t2, temp);
  fmpq_pow_si(temp, &ABCD[2], 2);
  fmpq_mul_si(temp, temp, 6);
  fmpq_add(t2, t2, temp);
  fmpq_mul(temp, &ABCD[0], &ABCD[3]);
  fmpq_mul_si(temp, temp, -3);
  fmpq_add(t2, t2, temp);

  fmpq_pow_si(t3, &ABCD[2], 2);
  fmpq_pow_si(temp, &ABCD[1], 3);
  fmpq_div_si(temp, temp, -6);
  fmpq_add(t3, t3, temp);

  res = (fmpq_is_zero(t1) && fmpq_is_zero(t2));
  if (res)
    {
      if (fmpq_is_zero(&ABCD[3]) || fmpq_is_zero(t3))
	{
	  flint_printf("(mestre_fmpz: bolza_conditions_23) Error: unexpected vanishing\n");
	  fflush(stdout);
	  flint_abort();
	  
	}
      /* Set t4 to denominator of alpha^2 */
      fmpq_pow_si(t4, &ABCD[0], 2);
      fmpq_mul(t4, t4, &ABCD[1]);
      fmpq_mul_si(t4, t4, 2);
      fmpq_mul(temp, &ABCD[0], &ABCD[2]);
      fmpq_mul_si(temp, temp, -3);
      fmpq_add(t4, t4, temp);
      fmpq_pow_si(temp, &ABCD[1], 2);
      fmpq_mul_si(temp, temp, -15);
      fmpq_add(t4, t4, temp);
      /* Set t1 to numerator */
      fmpq_pow_si(t1, &ABCD[1], 2);
      fmpq_mul(temp, &ABCD[0], &ABCD[2]);
      fmpq_add(t1, t1, temp);
      fmpq_div(t1, t1, t4);

      acb_set_fmpq(a, t1, prec);
      borchardt_sqrt(a, a, prec);
      acb_mul_si(a, a, 10, prec);
    }
  
  fmpq_clear(t1);
  fmpq_clear(t2);
  fmpq_clear(t3);
  fmpq_clear(t4);
  fmpq_clear(temp);
  return res;
}


int mestre_fmpz(acb_poly_t crv, fmpz* IC, slong prec)
{
  fmpq* ABCD;
  fmpq_t R2;
  acb_ptr IC_acb;
  acb_t a;
  slong j;
  int res = 1;

  ABCD = _fmpq_vec_init(4);
  fmpq_init(R2);
  IC_acb = _acb_vec_init(4);
  acb_init(a);
  for (j = 0; j < 4; j++) acb_set_fmpz(&IC_acb[j], &IC[j]);

  igusa_ABCD_from_IC_fmpz(ABCD, IC);
  igusa_R2_from_IC_fmpz(R2, IC);
  acb_poly_zero(crv);
  
  if (igusa_has_generic_automorphisms(IC_acb, prec))
    {
      res = mestre(crv, IC_acb, prec);
    }
  else if (fmpz_is_zero(&IC[3]))
    {
      flint_printf("(mestre_fmpz) Warning: not the invariants of a genus 2 curve\n");
      res = 0;
    }
  else if (!fmpq_is_zero(R2)) /* A,B,C are all zero: case II */
    {
      acb_poly_set_coeff_si(crv, 6, 1);
      acb_poly_set_coeff_si(crv, 1, 1);
    }
  /* R2 is zero */
  else if (bolza_conditions_12(ABCD)) /* Case VI */
    {
      acb_poly_set_coeff_si(crv, 5, 1);
      acb_poly_set_coeff_si(crv, 1, 1);
    }
  else if (bolza_conditions_11(ABCD))  /* Case V */
    {
      acb_poly_set_coeff_si(crv, 6, 1);
      acb_poly_set_coeff_si(crv, 0, 1);
    }
  else if (bolza_conditions_19(a, ABCD, prec)) /* Case IV */
    {
      acb_poly_set_coeff_acb(crv, 3, a);
      acb_poly_set_coeff_si(crv, 6, 1);
      acb_poly_set_coeff_si(crv, 0, 1);
    }
  else if (bolza_conditions_23(a, ABCD, prec)) /* Case III */
    {
      acb_poly_set_coeff_si(crv, 5, 1);
      acb_poly_set_coeff_acb(crv, 3, a);
      acb_poly_set_coeff_si(crv, 1, 1);
    }
  else /* Case I */
    {
      cardona(crv, IC_acb, prec);
    }

  _fmpq_vec_clear(ABCD, 4);
  fmpq_clear(R2);
  _acb_vec_clear(IC_acb, 4);
  acb_clear(a);
  return res;
}

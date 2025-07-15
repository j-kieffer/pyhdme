
#include "igusa.h"

static void igusa_ec_period_0(acb_t tau, slong prec)
{
  /* Return (1+isqrt(3))/2 */
  acb_one(tau);
  arb_sqrt_ui(acb_imagref(tau), 3, prec);
  acb_div_si(tau, tau, 2, prec);
}

static void igusa_ec_period_1728(acb_t tau, slong prec)
{
  /* Return i */
  acb_onei(tau);
}

static int igusa_ec_period_gen(acb_t tau, const acb_t j, slong prec)
{
  acb_t t, u;
  arb_t cmp;  
  acb_ptr kp2;
  acb_ptr k, kp;
  slong nb = 0;
  acb_mat_t m;
  int res = 1;
  slong i;

  acb_init(t);
  acb_init(u);
  arb_init(cmp);
  kp2 = _acb_vec_init(6);
  k = _acb_vec_init(6);
  kp = _acb_vec_init(6);
  acb_mat_init(m, 1, 1);

  res = igusa_ec_possible_kp2(kp2, j, prec);

  if (res)
    {
      for (i = 0; i < 6; i++)
	{
	  borchardt_sqrt(t, &kp2[i], prec);
	  acb_sub_si(&kp2[i], &kp2[i], 1, prec);
	  acb_neg(&kp2[i], &kp2[i]);
	  borchardt_sqrt(u, &kp2[i], prec);
      
	  /* Adjust signs given that real parts must be positive */
	  if (arb_is_negative(acb_realref(t))) acb_neg(t, t);
	  if (arb_is_negative(acb_realref(u))) acb_neg(u, u);
	  
	  /* If real part of k' contains zero, but not smaller than 0.1, 
	     we do not know what to do. */
	  arb_abs(cmp, acb_realref(t));
	  arb_mul_si(cmp, cmp, 10, prec);
	  arb_sub_si(cmp, cmp, 1, prec);
	  if (!arb_is_negative(cmp) && arb_contains_zero(acb_realref(t)))
	    {
	      res = 0;
	      break;
	    }
	  
	  /* If real part of k contains zero, but argument of k is not
	     outside [-1.5, 1.5], we do not know what to do. */
	  acb_arg(cmp, u, prec);
	  arb_abs(cmp, cmp);
	  arb_mul_si(cmp, cmp, 2, prec);
	  arb_sub_si(cmp, cmp, 3, prec);
	  if (!arb_is_positive(cmp) && arb_contains_zero(acb_realref(u)))
	    {
	      res = 0;
	      break;
	    }
	  
	  /* Now we can restrict to k, k' with positive real part */
	  if (arb_is_positive(acb_realref(t)) && arb_is_positive(acb_realref(u)))
	    {
	      acb_set(&kp[nb], t);
	      acb_set(&k[nb], u);
	      nb++;
	    }
	}
    }

  /* We have collected candidate values for k, k'. Now apply AGM,
     and stop when one of the results lies close to the fundamental domain */
  if (res)
    {
      res = 0;
      for (i = 0; i < nb; i++)
	{
	  acb_agm1(t, &kp[i], prec);
	  acb_agm1(u, &k[i], prec);
	  acb_div(t, t, u, prec);
	  acb_mul_onei(t, t);
	  acb_set(acb_mat_entry(m, 0, 0), t);
	  if (!siegel_not_in_fundamental_domain(m, prec))
	    {
	      res = 1;
	      acb_set(tau, t);
	      break;
	    }
	}
    }    

  acb_clear(t);
  acb_clear(u);
  arb_clear(cmp);
  _acb_vec_clear(kp2, 6);
  _acb_vec_clear(k, 6);
  _acb_vec_clear(kp, 6);
  acb_mat_clear(m); 
  return res;
}

int igusa_ec_period(acb_t tau, const acb_t j, slong prec)
{
  acb_t t;
  int res = 1;
  
  acb_init(t);
  acb_sub_si(t, j, 1728, prec);
  
  if (acb_equal_si(j, 0))
    {
      igusa_ec_period_0(tau, prec);
    }
  else if (acb_equal_si(j, 1728))
    {
      igusa_ec_period_1728(tau, prec);
    }
  else if (acb_contains_zero(j) || acb_contains_zero(t))
    {
      flint_printf("(igusa_ec_period) Warning: not sure if input is distinct from 0, 1728\n");
      res = 0;
    }
  else
    {
      res = igusa_ec_period_gen(tau, j, prec);
    }
  
  acb_clear(t);
  return res;
}


#include "theta.h"

/* prec is the common precision of tau and th_half. Target precision
   for derivatives is approximately prec/2. We approximate the
   derivatives to the value of finite differences.

   th_half is normalized with first entry 1. */

int theta_0123half_inverse_diff(acb_mat_t dtau, const acb_mat_t tau, acb_srcptr th_half,
				slong prec)
{
  acb_ptr th_half_pert;
  acb_mat_t tau_pert;
  arb_t eps;
  
  /* Try to maximize output precision */
  slong eps_exp = - (prec/2 - THETA_NEWTON_DERIVATIVE_OFFSET);
  slong g = 2;
  slong j;
  int res = 1;

  if (prec < THETA_NEWTON_MINPREC) return 0;

  acb_mat_init(tau_pert, g, g);
  th_half_pert = _acb_vec_init(4);
  arb_init(eps);

  arb_one(eps);
  arb_mul_2exp_si(eps, eps, eps_exp);

  for (j = 1; j < 4; j++)
    {
      /* Derivative wrt th_half[j] -> column j-1 of dtau */
      _acb_vec_set(th_half_pert, th_half, 4);
      acb_add_arb(&th_half_pert[j], &th_half_pert[j], eps, prec);
      res = res && theta_0123half_inverse_no_sqrt(tau_pert, th_half_pert, prec);
      /* Set tau_pert to contain [dtau_1, d(tau_3^2); d(tau_3^2), dtau_2] */
      acb_sub(acb_mat_entry(tau_pert, 0, 0), acb_mat_entry(tau_pert, 0, 0),
	      acb_mat_entry(tau, 0, 0), prec);
      acb_sub(acb_mat_entry(tau_pert, 1, 1), acb_mat_entry(tau_pert, 1, 1),
	      acb_mat_entry(tau, 1, 1), prec);
      
      acb_sqr(acb_mat_entry(tau_pert, 1, 0), acb_mat_entry(tau, 1, 0), prec);
      /* No need to square tau3 in tau_pert because we use
	 theta_0123half_inverse_no_sqrt */
					
      acb_sub(acb_mat_entry(tau_pert, 1, 0), acb_mat_entry(tau_pert, 0, 1),
	      acb_mat_entry(tau_pert, 1, 0), 2*prec);
      acb_set(acb_mat_entry(tau_pert, 0, 1), acb_mat_entry(tau_pert, 1, 0));
      
      acb_mat_scalar_mul_2exp_si(tau_pert, tau_pert, -eps_exp);
      acb_set(acb_mat_entry(dtau, 0, j-1), acb_mat_entry(tau_pert, 0, 0));
      acb_set(acb_mat_entry(dtau, 1, j-1), acb_mat_entry(tau_pert, 1, 1));
      acb_set(acb_mat_entry(dtau, 2, j-1), acb_mat_entry(tau_pert, 1, 0));
    }

  acb_mat_clear(tau_pert);
  _acb_vec_clear(th_half_pert, 4);
  arb_clear(eps);
  return res;
}

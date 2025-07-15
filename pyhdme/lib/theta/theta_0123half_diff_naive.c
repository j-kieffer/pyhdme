
#include "theta.h"

int theta_0123half_diff_naive(acb_mat_t dth, const acb_mat_t tau, slong prec)
{
  acb_mat_t tau_half;
  acb_ptr th_half;
  acb_ptr th_half_pert;
  arb_t eps;
  acb_t th0;
  mag_t err;
  
  /* Try to maximize output precision */
  slong eps_exp = - (prec/2 - THETA_NEWTON_DERIVATIVE_OFFSET);
  slong g = 2;
  slong i;
  int res = 1;

  acb_mat_init(tau_half, g, g);
  th_half = _acb_vec_init(4);
  th_half_pert = _acb_vec_init(4);
  arb_init(eps);
  acb_init(th0);
  mag_init(err);
  
  arb_one(eps);
  arb_mul_2exp_si(eps, eps, eps_exp);
  
  acb_mat_scalar_mul_2exp_si(tau_half, tau, -1);
  res = res && theta_0123_naive(th_half, tau_half, prec);
  acb_set(th0, &th_half[0]);
  _acb_vec_scalar_div(th_half, th_half, 4, th0, prec);

  /* First column of dth is derivative wrt tau1, second column wrt
     tau2, third column wrt tau3^2 */
  acb_mat_set(tau_half, tau);
  acb_add_arb(acb_mat_entry(tau_half, 0, 0), acb_mat_entry(tau_half, 0, 0), eps, prec);
  acb_mat_scalar_mul_2exp_si(tau_half, tau_half, -1);
  res = res && theta_0123_naive(th_half_pert, tau_half, prec);
  acb_set(th0, &th_half_pert[0]);
  _acb_vec_scalar_div(th_half_pert, th_half_pert, 4, th0, prec);
  _acb_vec_sub(th_half_pert, th_half_pert, th_half, 4, prec);
  _acb_vec_scalar_mul_2exp_si(th_half_pert, th_half_pert, 4, -eps_exp);
  for (i = 1; i < 4; i++)
    {
      acb_set(acb_mat_entry(dth, i-1, 0), &th_half_pert[i]);
    }
  
  acb_mat_set(tau_half, tau);
  acb_add_arb(acb_mat_entry(tau_half, 1, 1), acb_mat_entry(tau_half, 1, 1), eps, prec);
  acb_mat_scalar_mul_2exp_si(tau_half, tau_half, -1);
  res = res && theta_0123_naive(th_half_pert, tau_half, prec);
  acb_set(th0, &th_half_pert[0]);
  _acb_vec_scalar_div(th_half_pert, th_half_pert, 4, th0, prec);
  _acb_vec_sub(th_half_pert, th_half_pert, th_half, 4, prec);
  _acb_vec_scalar_mul_2exp_si(th_half_pert, th_half_pert, 4, -eps_exp);
  for (i = 1; i < 4; i++)
    {
      acb_set(acb_mat_entry(dth, i-1, 1), &th_half_pert[i]);
    }

  /* Set tau_3' to tau_3 \sqrt{1 + \eps/\tau_3^2} */
  acb_mat_set(tau_half, tau);
  acb_sqr(acb_mat_entry(tau_half, 0, 1), acb_mat_entry(tau_half, 0, 1), prec);
  acb_set_arb(acb_mat_entry(tau_half, 1, 0), eps);
  acb_div(acb_mat_entry(tau_half, 0, 1), acb_mat_entry(tau_half, 1, 0),
	  acb_mat_entry(tau_half, 0, 1), prec);
  acb_add_si(acb_mat_entry(tau_half, 0, 1), acb_mat_entry(tau_half, 0, 1), 1, prec);
  borchardt_sqrt(acb_mat_entry(tau_half, 0, 1), acb_mat_entry(tau_half, 0, 1), prec);
  acb_mul(acb_mat_entry(tau_half, 0, 1), acb_mat_entry(tau_half, 0, 1),
	  acb_mat_entry(tau, 0, 1), prec);
  acb_set(acb_mat_entry(tau_half, 1, 0), acb_mat_entry(tau_half, 0, 1));
  
  /* acb_add_arb(acb_mat_entry(tau_half, 0, 1), acb_mat_entry(tau_half, 0, 1), eps, prec);
     acb_add_arb(acb_mat_entry(tau_half, 1, 0), acb_mat_entry(tau_half, 1, 0), eps, prec); */
  acb_mat_scalar_mul_2exp_si(tau_half, tau_half, -1);
  res = res && theta_0123_naive(th_half_pert, tau_half, prec);
  acb_set(th0, &th_half_pert[0]);
  _acb_vec_scalar_div(th_half_pert, th_half_pert, 4, th0, prec);
  _acb_vec_sub(th_half_pert, th_half_pert, th_half, 4, prec);
  _acb_vec_scalar_mul_2exp_si(th_half_pert, th_half_pert, 4, -eps_exp);
  for (i = 1; i < 4; i++)
    {
      acb_set(acb_mat_entry(dth, i-1, 2), &th_half_pert[i]);
    }

  /* Add heuristic error */
  mag_set_ui_2exp_si(err, 1, eps_exp + THETA_DER_LOSS);
  acb_mat_add_error_mag(dth, err);
  
  acb_mat_clear(tau_half);
  _acb_vec_clear(th_half, 4);
  _acb_vec_clear(th_half_pert, 4);
  arb_clear(eps);
  acb_clear(th0);
  mag_clear(err);
  return res;
}

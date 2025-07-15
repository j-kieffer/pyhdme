
#include "theta.h"

int theta2_unif(acb_ptr th2, const acb_mat_t tau, slong prec)
{
  slong g = 2;
  slong k1 = 0;
  slong k2 = 0;

  fmpz_mat_t eta, eta_inv;
  fmpz_mat_t real_red;
  fmpz_mat_t M1;
  acb_ptr current_th2;
  acb_mat_t w;
  fmpz_t nb_real_red;
  arb_t tol;

  int flag_naive, flag_newton;
  int res;
  slong i, j;
  int v = get_theta_verbose();

  fmpz_mat_init(eta, 2*g, 2*g);
  fmpz_mat_init(eta_inv, 2*g, 2*g);
  fmpz_mat_init(real_red, 2*g, 2*g);
  fmpz_mat_init(M1, 2*g, 2*g);
  current_th2 = _acb_vec_init(n_pow(2, 2*g));
  acb_mat_init(w, g, g);
  fmpz_init(nb_real_red);
  arb_init(tol);

  arb_one(tol);
  arb_mul_2exp_si(tol, tol, THETA_NEWTON_TOL_EXP);

  res = siegel_is_in_fundamental_domain(tau, tol, prec);

  if (res)
    {
      k2 = theta_newton_k2(w, tau, prec);
      siegel_reduce_real(w, real_red, w, tol, prec);
      k1 = theta_newton_k1(w, w, prec);
      
      flag_naive = theta_use_naive(w, prec);
      flag_newton = theta_use_newton(w, prec);
      
      res = siegel_is_weakly_reduced(w, tol, prec)
	&& (flag_naive || flag_newton);
    }
    
  /* Compute theta2 at w */
  if (res && flag_naive)
    {
      if (v>0) {flint_printf("(theta2_unif) naive: %wd\n", prec); fflush(stdout);}
      res = theta2_naive(current_th2, w, prec);
    }
  else if (res && flag_newton)
    {
      if (v>0) {flint_printf("(theta2_unif) newton: %wd\n", prec); fflush(stdout);}
      res = theta2_newton(current_th2, w, prec);
    }

  /* Propagate along k1 sequence */
  for (i = 0; (i < k1) && res; i++)
    {
      /* Each step: rescale; compute square roots in four first
	 entries of current_th2; apply theta_duplication */
      for (j = 1; j < 4; j++)
	{
	  acb_div(&current_th2[j], &current_th2[j], &current_th2[0], prec);
	  res = res && acb_sqrt_goodpos(&current_th2[j], &current_th2[j], prec);
	}
      acb_one(&current_th2[0]);
      /* theta_duplication uses only the first 4 entries and allows
	 aliasing */
      theta_duplication(current_th2, current_th2, prec);
    }
  
  /* We're back to the end of the k2 sequence. Propagate again */
  fmpz_mat_zero(eta);
  fmpz_set_si(fmpz_mat_entry(eta, 0, 2), 1);
  fmpz_set_si(fmpz_mat_entry(eta, 1, 1), 1);
  fmpz_set_si(fmpz_mat_entry(eta, 2, 0), -1);
  fmpz_set_si(fmpz_mat_entry(eta, 3, 3), 1);
  fmpz_mat_direct_inv(eta_inv, eta);
  fmpz_mat_one(M1);
  
  /* Real reduction happens only on the y1 coordinate */
  fmpz_set(nb_real_red, fmpz_mat_entry(real_red, 0, 2));
  if (arb_is_positive(acb_realref(acb_mat_entry(tau, 0, 0))))
    {
      fmpz_set_si(fmpz_mat_entry(M1, 0, 2), 1);
    }
  else
    {
      fmpz_set_si(fmpz_mat_entry(M1, 0, 2), -1);
    }
  
  for (i = 0; (i < k2) && res; i++)
    {
      /* Each step: apply translation if nb_real_red is currently odd;
       apply theta_transform with eta; square roots and duplication;
       finally apply theta2_transform with eta_inv. */
      if (fmpz_is_odd(nb_real_red))
	{
	  theta2_transform(current_th2, M1, current_th2, prec);
	}
      /* Divide by 2, round towards zero */
      fmpz_tdiv_q_2exp(nb_real_red, nb_real_red, 1);
      /* Get theta constants at eta * current */
      theta2_transform(current_th2, eta, current_th2, prec);
      /* Rescale and square roots */
      for (j = 1; j < 4; j++)
	{
	  acb_div(&current_th2[j], &current_th2[j], &current_th2[0], prec);
	  res = res && acb_sqrt_goodpos(&current_th2[j], &current_th2[j], prec);
	}
      acb_one(&current_th2[0]);
      /* Duplication */
      theta_duplication(current_th2, current_th2, prec);
      /* Apply eta_inv */
      theta2_transform(current_th2, eta_inv, current_th2, prec);
    }

  /* Store end result */
  _acb_vec_set(th2, current_th2, n_pow(2, 2*g));

  /* Print warning */
  if (!res && (v>0))
    {
      flint_printf("(theta2_unif) Warning: computation failed for the matrix\n");
      acb_mat_printd(tau, 10);
    }
  
  fmpz_mat_clear(eta);
  fmpz_mat_clear(eta_inv);
  fmpz_mat_clear(real_red);
  fmpz_mat_clear(M1);
  _acb_vec_clear(current_th2, n_pow(2, 2*g));
  acb_mat_clear(w);
  fmpz_clear(nb_real_red);
  arb_clear(tol);

  return res;
}

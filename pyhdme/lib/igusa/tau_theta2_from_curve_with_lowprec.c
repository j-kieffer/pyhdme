
#include "igusa.h"

/* See also tau_theta2_from_curve */

int tau_theta2_from_curve_with_lowprec(acb_mat_t tau, acb_ptr theta2, const acb_poly_t crv,
    acb_srcptr th2_lp, slong prec)
{
  time_pair start; timestamp_mark(&start);
  slong perm, signs;
  acb_ptr I;
  acb_ptr roots;
  acb_ptr ros;
  acb_ptr th2, th4;
  int res;
  int v = get_thomae_verbose();

  I = _acb_vec_init(4);
  roots = _acb_vec_init(6);
  ros = _acb_vec_init(3);
  th2 = _acb_vec_init(16);
  th4 = _acb_vec_init(16);

  igusa_from_curve(I, crv, prec);

  if (v) flint_printf("(tau_theta2_from_curve_with_lowprec) Looking for correct root ordering...\n");
  res = thomae_roots(roots, crv, prec);

  if (res) res = thomae_correct_signs_with_lowprec(&perm, &signs, roots, I, th2_lp, prec);
  if (res) {
    time_pair mid; timestamp_mark(&mid);
    if(v) flint_printf("(tau_theta2_from_curve_with_lowprec:from roots) Computing period matrix at high precision from roots...\n");
    thomae_reorder(roots, roots, perm);
    thomae_rosenhain(ros, roots, prec);
    thomae_theta4(th4, ros, prec);
    thomae_theta2(th2, th4, ros, signs, prec);
    res = theta2_inverse(tau, th2, prec);
    if (v) report("tau_theta2_from_curve_with_lowprec:from roots", mid);
  }

  theta2_renormalize(theta2, th2, prec);

  _acb_vec_clear(I, 4);
  _acb_vec_clear(roots, 6);
  _acb_vec_clear(ros, 3);
  _acb_vec_clear(th2, 16);
  _acb_vec_clear(th4, 16);
  if (v) report_end(start);
  return res;
}

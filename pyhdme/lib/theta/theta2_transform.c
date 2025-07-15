
#include "theta.h"

/* NB: In the context of the theta2 function, we should perhaps
   compute once and for all the transformations associated with the
   two eta matrices we use. */

void theta2_transform(acb_ptr th2_eta, const fmpz_mat_t eta, acb_srcptr th2, slong prec)
{
  acb_ptr res;
  slong g = 2;
  ulong n = n_pow(2, 2*g);
  ulong ch;
  ulong image_ch;
  fmpz_t epsilon;
  acb_t th2_eta_entry;

  res = _acb_vec_init(n);
  fmpz_init(epsilon);
  acb_init(th2_eta_entry);

  for (ch = 0; ch < n; ch++)
    {
      image_ch = theta_transform_image_char(epsilon, ch, eta);
      acb_unit_root(th2_eta_entry, 4, prec); /* Only change wrt theta_transform */
      acb_pow_fmpz(th2_eta_entry, th2_eta_entry, epsilon, prec);
      acb_mul(th2_eta_entry, th2_eta_entry, &th2[image_ch], prec);
      acb_set(&res[ch], th2_eta_entry);
    }

  _acb_vec_set(th2_eta, res, n);

  _acb_vec_clear(res, n);
  fmpz_clear(epsilon);
  acb_clear(th2_eta_entry);
}

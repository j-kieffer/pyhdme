
#include "theta.h"

/* NB: In the context of the theta2 function, we should perhaps
   compute once and for all the transformations associated with the
   two eta matrices we use. */

void theta_transform(acb_ptr th_eta, const fmpz_mat_t eta, acb_srcptr th, slong prec)
{
  acb_ptr res;
  slong g = 2;
  ulong n = n_pow(2, 2*g);
  ulong ch;
  ulong image_ch;
  fmpz_t epsilon;
  acb_t th_eta_entry;

  res = _acb_vec_init(n);
  fmpz_init(epsilon);
  acb_init(th_eta_entry);

  for (ch = 0; ch < n; ch++)
    {
      image_ch = theta_transform_image_char(epsilon, ch, eta);
      acb_unit_root(th_eta_entry, 8, prec);
      acb_pow_fmpz(th_eta_entry, th_eta_entry, epsilon, prec);
      acb_mul(th_eta_entry, th_eta_entry, &th[image_ch], prec);
      acb_set(&res[ch], th_eta_entry);
    }

  _acb_vec_set(th_eta, res, n);

  _acb_vec_clear(res, n);
  fmpz_clear(epsilon);
  acb_clear(th_eta_entry);
}

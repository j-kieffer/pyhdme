
#include "theta.h"

void theta_transform_matrix(fmpz_mat_t res, const fmpz_mat_t eta)
{
  slong label;
  ulong ch;
  slong image_label;
  ulong image_ch;
  fmpz_t epsilon;
  slong g = 2;

  fmpz_init(epsilon);

  fmpz_mat_zero(res);

  for (ch = 0; ch < n_pow(2, 2*g); ch++)
    {
      if (theta_char_is_even(ch, g))
	{
	  label = theta_char_get_label_g2(ch);
	  image_ch = theta_transform_image_char(epsilon, ch, eta);
	  image_label = theta_char_get_label_g2(image_ch);

	  fmpz_set_si(fmpz_mat_entry(res, label, 0), label);
	  fmpz_set(fmpz_mat_entry(res, label, 1), epsilon);
	  fmpz_set_si(fmpz_mat_entry(res, label, 2), image_label);
	}
    }
  
  fmpz_clear(epsilon);
}

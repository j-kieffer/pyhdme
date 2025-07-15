
#include "modular.h"

int modeq_round(modeq_t R, slong* bits, const modeq_acb_t E)
{  
  arf_t radius, max_radius;
  slong radius_bits;
  int b;
  int res = 1;
  slong k;
  slong v = get_modeq_verbose();

  arf_init(radius);
  arf_init(max_radius);

  modeq_degree(R) = modeq_degree(E);
  modeq_nb(R) = modeq_nb(E);

  res = acb_round(modeq_den(R), max_radius, modeq_den(E));
  
  for (k = 0; k < modeq_nb(E)+1; k++)
    {
      b = acb_poly_round(&modeq_all_nums(R)[k], radius,
			 &modeq_all_nums(E)[k], modeq_degree(E));
      res = res && b;
      arf_max(max_radius, max_radius, radius);
    }

  radius_bits = arf_abs_bound_lt_2exp_si(max_radius);
  *bits = radius_bits;
  
  if (v && res) flint_printf("(modeq_round) Excess precision: %wd\n", - *bits);
  if (v && !res) flint_printf("(modeq_round) Need %wd more bits of precision\n", *bits);
  if (res) pol_simplify(modeq_all_nums(R), modeq_den(R), modeq_degree(R),
			modeq_nb(R)+1);
  
  arf_clear(radius);
  arf_clear(max_radius);
  return res;  
}

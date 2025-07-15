
#include "igusa.h"

int igusa_has_generic_automorphisms(acb_srcptr IC, slong prec)
{
  acb_ptr ABCD;
  acb_t R2;
  int res;

  ABCD = _acb_vec_init(4);
  acb_init(R2);

  igusa_R2_from_IC(R2, IC, prec);
  igusa_ABCD_from_IC(ABCD, IC, prec);

  res = !acb_contains_zero(&IC[3])
    && !acb_contains_zero(R2)
    && (!acb_contains_zero(&ABCD[0])
	|| !acb_contains_zero(&ABCD[1])
	|| !acb_contains_zero(&ABCD[2]));
  /*
  if (!res)
    {
      flint_printf("(igusa_has_generic_automorphisms) Warning: cannot exclude extra automorphisms\n");
      flint_printf("(igusa_has_generic_automorphisms) R2: "); acb_printd(R2, 10);
      flint_printf("\n(igusa_has_generic_automorphisms) Clebsch:\n");
      acb_printd(&ABCD[0], 10); flint_printf("\n");
      acb_printd(&ABCD[1], 10); flint_printf("\n");
      acb_printd(&ABCD[2], 10); flint_printf("\n");
      } */

  _acb_vec_clear(ABCD, 4);
  acb_clear(R2);
  return res;
}


#include "theta.h"

ulong theta_char_set_label_g2(slong label)
{
  /* label is a1 a0 b1 b0, we want a0 a1 b0 b1 */
  ulong a = n_revbin(label>>2, 2);
  ulong b = n_revbin(label%4, 2);
  return theta_char_set_ab(a, b, 2);
}


#include "theta.h"

slong theta_char_get_label_g2(ulong ch)
{
  ulong a = theta_char_get_a(ch, 2);
  ulong b = theta_char_get_b(ch, 2);
  return (n_revbin(a, 2)<<2) + n_revbin(b, 2);
}

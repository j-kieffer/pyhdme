
#include "theta.h"

int
theta_char_is_even(ulong ch, slong g)
{
  ulong a = theta_char_get_a(ch, g);
  ulong b = theta_char_get_b(ch, g);
  int res = (theta_char_dot_product(a, b, g) % 2 == 0);
  return res;
}

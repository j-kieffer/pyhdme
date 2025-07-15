
#include "theta.h"

/* Assumes that the ch are all distinct. */

int theta_char_is_goepel(ulong ch1, ulong ch2, ulong ch3, ulong ch4, slong g)
{
  return theta_char_is_even(ch1, g)
    && theta_char_is_even(ch2, g)
    && theta_char_is_even(ch3, g)
    && theta_char_is_even(ch4, g)
    && ((ch1 ^ ch2 ^ ch3 ^ ch4) == 0);
}


#include "theta.h"

/* Assumes that the ch are all distinct. */

/* That's inefficient: compute ch from the Goepel relation, then check
   if it is even, etc. */

int theta_char_is_syzygous(ulong ch1, ulong ch2, ulong ch3, slong g)
{
  ulong n = n_pow(2, 2*g);
  int res = 0;
  ulong ch;

  for (ch = 0; ch < n; ch++)
    {
      if (theta_char_is_even(ch, g)
	  && (ch != ch1)
	  && (ch != ch2)
	  && (ch != ch3)
	  && theta_char_is_goepel(ch, ch1, ch2, ch3, g))
	{
	  res = 1;
	}
    }
  return res;
}

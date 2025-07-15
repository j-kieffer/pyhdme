
#include "theta.h"

int
theta_char_dot_product(ulong a, ulong b, slong g)
{
  int res = 0;
  ulong and = a & b;
  while (and > 0)
    {
      if ((and & 1) == 1) res++;
      and = and>>1;
    }
  return res;
}



#include "theta.h"

ulong
theta_char_get_b(ulong ch, slong g)
{
  ulong mask = (1<<g) - 1;
  return ch & mask;
}

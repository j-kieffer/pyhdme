
#include "igusa.h"

slong thomae_startprec(slong prec)
{
  slong res = prec;
  while (res > THOMAE_MULPREC * (THOMAE_LOWPREC-1))
    {
      res = (res/THOMAE_MULPREC) + 1;
    }
  return res;
}

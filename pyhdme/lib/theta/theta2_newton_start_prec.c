
#include "theta.h"

slong theta2_newton_start_prec(slong target_prec)
{
  slong res = target_prec;
  while (res > THETA_NEWTON_BASEPREC)
    {
      res += THETA_NEWTON_LOSS;
      res = (res+1)/2;
    }
  return res;
}

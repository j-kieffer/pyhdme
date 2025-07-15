
#include "theta.h"

void theta_der_set_pert(arb_t eps, slong prec)
{
  slong eps_exp = - (prec/2 - THETA_NEWTON_DERIVATIVE_OFFSET);
  arb_one(eps);
  arb_mul_2exp_si(eps, eps, eps_exp);
}

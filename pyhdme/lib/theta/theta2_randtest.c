
#include "theta.h"

void theta2_randtest(acb_ptr theta2, flint_rand_t state, slong prec)
{
  acb_ptr th_half;

  th_half = _acb_vec_init(4);
  acb_one(&th_half[0]);
  acb_randtest_precise(&th_half[1], state, prec, 1);
  acb_randtest_precise(&th_half[2], state, prec, 1);
  acb_randtest_precise(&th_half[3], state, prec, 1);
  theta_duplication(theta2, th_half, prec);

  _acb_vec_clear(th_half, 4);
}

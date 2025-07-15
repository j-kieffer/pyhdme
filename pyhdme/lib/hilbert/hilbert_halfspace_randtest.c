
#include "hilbert.h"

void hilbert_halfspace_randtest(acb_ptr t, flint_rand_t state, slong prec)
{
  slong mag_bits = 1 + n_randint(state, 5);
  int stop;
  slong k;

  for (k = 0; k < 2; k++)
    {
      arb_randtest_precise(acb_realref(&t[k]), state, prec, mag_bits);
      stop = 0;
      while (!stop)
	{
	  arb_randtest_precise(acb_imagref(&t[k]), state, prec, mag_bits);
	  stop = arb_is_positive(acb_imagref(&t[k]));
	}
    }
}

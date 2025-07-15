
#include "modular.h"

void modeq_rescale(modeq_acb_t R, const modeq_acb_t E,
		   const acb_t c, slong prec)
{
  slong j;
  
  modeq_acb_set(R, E);

  acb_poly_scalar_mul(modeq_equation(R), modeq_equation(R), c, prec);
  for (j = 0; j < modeq_nb(R); j++)
    {
      acb_poly_scalar_mul(modeq_interpolate(R, j), modeq_interpolate(R, j),
			  c, prec);
    }
  acb_mul(modeq_den(R), modeq_den(R), c, prec);
}

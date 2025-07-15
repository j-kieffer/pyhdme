
#include "modular.h"

void modeq_acb_set(modeq_acb_t R, const modeq_acb_t E)
{
  slong j;
  
  modeq_degree(R) = modeq_degree(E);
  modeq_nb(R) = modeq_nb(E);
  
  acb_poly_set(modeq_equation(R), modeq_equation(E));
  for (j = 0; j < modeq_nb(E); j++)
    {
      acb_poly_set(modeq_interpolate(R, j), modeq_interpolate(E, j));
    }
  acb_set(modeq_den(R), modeq_den(E));  
}

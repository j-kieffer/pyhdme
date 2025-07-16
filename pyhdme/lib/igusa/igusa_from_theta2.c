
#include "igusa.h"

void igusa_from_theta2(acb_ptr I, acb_srcptr theta2, slong prec)
{
    acb_theta_g2_even_weight(&I[0], &I[1], &I[2], &I[3], theta2, prec);
}

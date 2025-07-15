
#include "modular.h"

slong siegel_direct_isog_startprec(fmpz* I, slong ell)
{
    slong weights[4] = IGUSA_WEIGHTS;
    slong h = cov_height(I, 4, weights);

    slong res = 10*h + 30*n_clog(ell,2) + 90;

    return 100 * (res/100 + 1);

}

/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#include <flint/acb_mat.h>
#include <flint/acb_theta.h>
#include "theta.h"

int theta2_unif(acb_ptr th2, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_ptr z;

    z = _acb_vec_init(g);
    acb_theta_all(th2, z, tau, 1, prec);
    _acb_vec_clear(z, g);
    return 1;
}

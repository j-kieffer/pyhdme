/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef SIEGEL_H
#define SIEGEL_H

#include <flint/flint.h>
#include <flint/arb_mat.h>
#include <flint/acb_mat.h>
#include "flint_compat.h"

int siegel_not_minkowski_reduced(const arb_mat_t r, slong prec);

int siegel_not_real_reduced(const acb_mat_t z, slong prec);

int siegel_not_in_fundamental_domain(const acb_mat_t z, slong prec);

#endif

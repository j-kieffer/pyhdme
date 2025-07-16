/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef THETA_H
#define THETA_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/ulong_extras.h>

#include "flint_compat.h"
#include "siegel.h"
#include "verbose.h"

/* Borchardt means */

int acb_sqrt_goodpos(acb_t r, const acb_t z, slong prec);

void borchardt_sqrt(acb_t r, const acb_t z, slong prec);

void borchardt_root_ui(acb_t r, const acb_t z, ulong e, slong prec);

int borchardt_step(acb_ptr b, acb_srcptr a, slong prec);

void borchardt_mean_m0(arb_t m0, acb_srcptr a, slong prec);

void borchardt_mean_M0(arb_t M0, acb_srcptr a, slong prec);

void borchardt_mean_Delta0(arb_t Delta0, acb_srcptr a, slong prec);

int borchardt_mean_nb_steps_before_quad_conv(fmpz_t nb, acb_srcptr a, slong prec);

void borchardt_mean_nb_steps_after_quad_conv(fmpz_t nb, acb_srcptr a, slong prec);

int borchardt_mean_quad_conv_is_reached(acb_srcptr a, slong prec);

int borchardt_mean(acb_t r, acb_srcptr a, slong prec);

void borchardt_excl_half_planes(arf_struct* b, const acb_t z, slong prec);

int borchardt_mean_invalid(acb_srcptr a, slong prec);


/* Theta characteristics */

ulong theta_char_get_a(ulong ch, slong g);

ulong theta_char_get_b(ulong ch, slong g);

ulong theta_char_set_ab(ulong a, ulong b, slong g);

int theta_char_dot_product(ulong a, ulong b, slong g);

int theta_char_is_even(ulong ch, slong g);

slong theta_char_get_label_g2(ulong ch);

ulong theta_char_set_label_g2(slong label);


/* Theta constants */

void theta_duplication(acb_ptr th2_2tau, acb_srcptr th_tau, slong prec);

int theta2_inverse(acb_mat_t tau, acb_srcptr th, slong prec);

int theta2_invalid(acb_srcptr th2, slong prec);

int theta2_unif(acb_ptr th2, const acb_mat_t tau, slong prec);

int theta2_renormalize(acb_ptr th2, acb_srcptr th2_proj, slong prec);

void theta2_randtest(acb_ptr theta2, flint_rand_t state, slong prec);


#endif

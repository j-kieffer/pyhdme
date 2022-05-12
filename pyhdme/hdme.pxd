# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.

from sage.libs.flint.types cimport fmpq, slong

cdef extern from "slong.h":
    pass

cdef extern from "lib/modular.h":
    int siegel_modeq_isog_invariants_Q(
        slong* nb_roots,
        fmpq* all_isog_j,
        fmpq* j,
        slong ell)

    int siegel_modeq_2step_isog_invariants_Q(
        slong* nb_roots,
        fmpq* all_isog_j,
        fmpq* j,
        slong ell)


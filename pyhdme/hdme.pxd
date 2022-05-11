# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.

from sage.libs.flint.types cimport fmpq, mp_limb_signed_t


cdef extern from "lib/modular.h":
    int siegel_modeq_isog_invariants_Q(
        mp_limb_signed_t* nb_roots,
        fmpq* all_isog_j,
        fmpq* j,
        mp_limb_signed_t ell)

#    int siegel_modeq_2step_isog_invariants_Q(
#        slong* nb_roots,
#        fmpq* all_isog_j,
#        fmpq* j,
#        slong ell)

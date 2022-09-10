# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.

from sage.libs.flint.types cimport fmpz, slong

cdef extern from "slong.h":
    pass


cdef extern from "lib/modular.h":
    int siegel_direct_isog_Q(slong* nb, fmpz* all_I, fmpz* I, slong ell);
    int siegel_2step_direct_isog_Q(slong* nb, fmpz* all_I, fmpz* I, slong ell);

cdef extern from "lib/igusa.h":
    void igusa_IC_fmpz(fmpz* IC, fmpz* I);
    void igusa_from_IC_fmpz(fmpz* I, fmpz* IC);


cdef extern from "lib/verbose.h":
    int set_modeq_verbose(int i)
    int set_hecke_verbose(int i)
    int set_thomae_verbose(int i)



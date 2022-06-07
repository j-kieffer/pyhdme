# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.

from sage.libs.flint.types cimport slong, fmpz_poly_t

cdef extern from "slong.h":
    pass

cdef extern from "lib/hecke.h":
    int hecke_charpoly(
    	fmpz_poly_t pol,
	slong ell,
	slong wt)


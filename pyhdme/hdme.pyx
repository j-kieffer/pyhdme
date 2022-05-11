# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.


from cysignals.signals cimport sig_on, sig_off
from sage.libs.flint.types cimport fmpq
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational


def something_else(
    j,
    ell
):
    cdef slong* nb_roots = calloc(sizeof(slong))
    cdef fmpq* all_isog_j = calloc((ell^3 + ell^2 + ell + 1)*sizeof(fmpq))
    cdef fmpq* cj = calloc(sizeof(fmpq))
    cdef slong cell = slong(ell)
    for i in range(ell^3 + ell^2 + ell + 1):
        fmpq_init(all_isog_j[i])
    fmpq_init(cj)
    fmpq_set_mpq(cj, j.value)

    sig_on()
    siegel_modeq_isog_invariants_Q(nb_roots,
                                   all_isog_j,
                                   j,
                                   cell)
    sig_off()

	# convert all_isog_j to list of rationals:
	cdef Rational read = <Rational>Rational.__new__(Rational)
	cdef Integer nb_roots_python = Integer(nb_roots)
	ans = [Rational(0) for i in range(3*nb_roots_python)]
	for i in range(3*nb_roots_python):
		fmpt_get_mpq(read.value, (*all_isog_j)[i])
		ans[i] = read
	for i in range(ell^3 + ell^2  +ell + 1):
        fmpq_clear(all_isog_j[i])
	fmpq_clear(cj)
	cdef free(nb_roots)
	cdef free(all_isog_j)
	cdef free(cj)
    return ans


def siegel_modeq_2step_isog_invariants_Q():
    pass

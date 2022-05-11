# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.


from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator
from sage.libs.flint.types cimport fmpq
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational


def something_else(
    j,
    ell
):
    assert len(j) == 3
    assert ell.is_prime()
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef slong nb_roots;
    cdef fmpq* all_isog_j = calloc(ell^3 + ell^2 + ell + 1, sizeof(fmpq))
    cdef fmpq* cj = calloc(3, sizeof(fmpq))
    cdef slong cell = slong(ell)
    for i in range(ell^3 + ell^2 + ell + 1):
        fmpq_init(all_isog_j[i])
    fmpq_init(cj)
    for i in range(3):
        fmpq_set_mpq(cj + i, j[i].value)

    sig_on()
    siegel_modeq_isog_invariants_Q(nb_roots, all_isog_j, j, cell)
    sig_off()

	# convert all_isog_j to list of rationals:
	nb_roots_python = nb_roots.value
	ans = [Rational(0) for i in range(3*nb_roots_python)]
	for i in range(3*nb_roots_python):
		fmpt_get_mpq((<Rational>ans[i]).value, (*all_isog_j)[i])
		ans[i] = read
	for i in range(ell^3 + ell^2  +ell + 1):
        fmpq_clear(all_isog_j[i])
	fmpq_clear(cj)
    return ans


def siegel_modeq_2step_isog_invariants_Q():
    pass

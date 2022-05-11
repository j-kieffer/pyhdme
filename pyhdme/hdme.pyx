# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.


from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator
from sage.libs.flint.types cimport fmpq
from sage.libs.flint.fmpq cimport fmpq_init, fmpq_clear, fmpq_get_mpq, fmpq_set_mpq
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from cpython.int cimport PyInt_FromLong


def something_else(
    j,
    ell
):
    assert len(j) == 3
    assert ell.is_prime()
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef slong nb_roots;
    cdef fmpq* all_isog_j = mem.calloc(ell^3 + ell^2 + ell + 1, sizeof(fmpq))
    cdef fmpq* cj = mem.calloc(3, sizeof(fmpq))
    cdef slong cell = slong(ell)
    for i in range(3):
        fmpq_init(cj + i)
        fmpq_set_mpq(cj + i, j[i].value)

    sig_on()
    siegel_modeq_isog_invariants_Q(nb_roots, all_isog_j, j, cell)
    sig_off()

    # convert all_isog_j to list of rationals:
    cdef object nb_roots_python
    nb_roots_python = PyInt_FromLong(nb_roots)
    ans = [Rational(0) for _ in range(3*nb_roots_python)]
    for i in range(3 * nb_roots_python):
        fmpq_get_mpq((<Rational>ans[i]).value, &all_isog_j[i])

    for i in range(3):
        fmpq_clear(&cj[i])
    return ans


def siegel_modeq_2step_isog_invariants_Q():
    pass
# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.


from sage.libs.flint.types cimport fmpq, slong
from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator
from sage.libs.flint.fmpq cimport fmpq_init, fmpq_clear, fmpq_get_mpq, fmpq_set_mpq
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from cpython.int cimport PyInt_FromLong
from sage.libs.gmp.mpz cimport mpz_get_si
from sage.all import QQ, ZZ


def siegel_modeq_isog_invariants_Q_wrapper(
    absolute_invariants,
    ell,
    verbose=False
):
    # cdef int cverbose = int(verbose)
    # set_modeq_verbose(cverbose);
    assert len(absolute_invariants) == 3
    j = [QQ(elt) for elt in absolute_invariants]
    ell = ZZ(ell)
    assert ell.is_prime()
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef slong nb_roots;
    max_nb_roots = ell**3 + ell**2 + ell + 1
    cdef fmpq* all_isog_j = <fmpq*>mem.calloc(3 * max_nb_roots, sizeof(fmpq))
    cdef fmpq* cj = <fmpq*>mem.calloc(3, sizeof(fmpq))
    cdef slong cell = mpz_get_si((<Integer>ell).value)
    for i in range(3):
        fmpq_init(&cj[i])
        fmpq_set_mpq(&cj[i], (<Rational>j[i]).value)
    #for i in range(3*max_nb_roots):
    #    fmpq_init(&all_isog_j[i]);

    sig_on()
    assert siegel_modeq_isog_invariants_Q(&nb_roots, all_isog_j, cj, cell) == 1
    sig_off()

    # convert all_isog_j to list of rationals:
    cdef object nb_roots_python
    nb_roots_python = PyInt_FromLong(nb_roots)
    ans = [Rational(0) for _ in range(3*nb_roots_python)]
    for i in range(3 * nb_roots_python):
        fmpq_get_mpq((<Rational>ans[i]).value, &all_isog_j[i])

    for i in range(3):
        fmpq_clear(&cj[i])
    #for i in range(3*max_nb_roots):
    #    fmpq_clear(&all_isog_j[i]);
    res = [(ans[3*i], ans[3*i + 1], ans[3*i + 2]) for i in range(nb_roots_python)]
    print(res)
    return res


def siegel_modeq_2step_isog_invariants_Q():
    pass

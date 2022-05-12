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



def igusa_clebsch_to_hdme_absolute_invariants(I2, I4, I6, I10):
    I6prime = (I2*I4-3*I6)/2
    return (I4*I6prime/I10, I2*I4*I4/I10, I4**5/I10**2)


def hdme_absolute_invariants_to_igusa_clebsch(j1, j2, j3):
    # M = Matrix([[0,1,1,-1],[1,2,0,-1],[0,5,0,-2]])
    # M.solve_left(Matrix([[-2,1,0,0],[-3,0,1,0],[-5,0,0,1]]))
    # [ 0 -2  1]
    # [ 1 -3  1]
    # [ 0 -5  2]
    I2, I4, I6prime, I10 = (1, j3/j2**2, j1*j3/j2**3, j3**2/j2**5)
    I6 = (I2 * I4 - 2*I6prime)/3
    return (I2, I4, I6, I10)


def generic_wrapper(
    steps,
    absolute_invariants,
    ell,
    verbose=False
):
    # cdef int cverbose = int(verbose)
    # set_modeq_verbose(cverbose);
    assert steps in [1,2]
    assert len(absolute_invariants) == 3
    j = [QQ(elt) for elt in absolute_invariants]
    ell = ZZ(ell)
    assert ell.is_prime()
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef slong nb_roots;
    max_nb_roots = (ell**3 + ell**2 + ell + 1)**steps
    cdef fmpq* all_isog_j = <fmpq*>mem.calloc(3 * max_nb_roots, sizeof(fmpq))
    cdef fmpq* cj = <fmpq*>mem.calloc(3, sizeof(fmpq))
    cdef slong cell = mpz_get_si((<Integer>ell).value)
    for i in range(3):
        fmpq_init(&cj[i])
        fmpq_set_mpq(&cj[i], (<Rational>j[i]).value)
    for i in range(3*max_nb_roots):
        fmpq_init(&all_isog_j[i]);

    sig_on()
    if steps == 1:
        assert siegel_modeq_isog_invariants_Q(&nb_roots, all_isog_j, cj, cell) == 1
    elif steps == 2:
        assert siegel_modeq_2step_isog_invariants_Q(&nb_roots, all_isog_j, cj, cell) == 1
    sig_off()

    # convert all_isog_j to list of rationals:
    cdef object nb_roots_python
    nb_roots_python = PyInt_FromLong(nb_roots)
    assert nb_roots_python <= max_nb_roots;
    ans = [Rational(0) for _ in range(3*nb_roots_python)]
    for i in range(3 * nb_roots_python):
        fmpq_get_mpq((<Rational>ans[i]).value, &all_isog_j[i])

    for i in range(3):
        fmpq_clear(&cj[i])
    for i in range(3*max_nb_roots):
        fmpq_clear(&all_isog_j[i]);
    res = [(ans[3*i], ans[3*i + 1], ans[3*i + 2]) for i in range(nb_roots_python)]
    return res


def siegel_modeq_2step_isog_invariants_Q_wrapper(
    igusa_clebsch_invariants,
    ell,
    verbose=False
):
    r"""
    Given the Igusa-Clebsch invariants `I_2, I_4, I_6, I_{10}` of Igusa and Clebsch [IJ1960]_
    of an abelian surface compute the Igusa-Clebsch invariants of surfaces (ell,ell)^2-isogenous

    INPUT:

    - ``igusa_clebsch_invariants`` -- the Igusa-Clebsch invariants `I_2, I_4, I_6, I_{10}` of Igusa and Clebsch [IJ1960]_
    - ``ell`` -- a prime number

    OUTPUT:

    - A list of Igusa-Clebsch invariants of abelian surfaces (ell, ell)-isogenous

    Examples::


    """
    I2, I4, I6, I10 = igusa_clebsch_invariants
    absolute_invariants = igusa_clebsch_to_hdme_absolute_invariants(I2, I4, I6, I10)
    res = generic_wrapper(1, absolute_invariants, ell, verbose=verbose)
    return [hdme_absolute_invariants_to_igusa_clebsch(*elt) for elt in res]

def siegel_modeq_isog_invariants_Q_wrapper(
    igusa_clebsch_invariants,
    ell,
    verbose=False
):
    r"""
    Given the Igusa-Clebsch invariants `I_2, I_4, I_6, I_{10}` of Igusa and Clebsch [IJ1960]_
    of an abelian surface compute the Igusa-Clebsch invariants of surfaces (ell,ell)-isogenous

    INPUT:

    - ``igusa_clebsch_invariants`` -- the Igusa-Clebsch invariants `I_2, I_4, I_6, I_{10}` of Igusa and Clebsch [IJ1960]_
    - ``ell`` -- a prime number

    OUTPUT:

    - A list of Igusa-Clebsch invariants of abelian surfaces (ell, ell)-isogenous

    Examples::


    """
    I2, I4, I6, I10 = igusa_clebsch_invariants
    absolute_invariants = igusa_clebsch_to_hdme_absolute_invariants(I2, I4, I6, I10)
    res = generic_wrapper(1, absolute_invariants, ell, verbose=verbose)
    return [hdme_absolute_invariants_to_igusa_clebsch(*elt) for elt in res]

